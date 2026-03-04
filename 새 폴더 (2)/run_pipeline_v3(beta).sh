#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════
# MTG2 한우 4형질 다변량 유전체 육종가 추정 파이프라인 (v2)
# ═══════════════════════════════════════════════════════════════════════
#
# [특징]
#   1. ID: GRM의 원래 ID를 그대로 사용 (리넘버링 없음)
#   2. -cc (분류형 공변량), -qc (연속형 공변량) 옵션 지원
#   3. 4형질 동시 분석 (다변량 모델, -mod 4 -cove 1)
#
# [필수 입력 파일]
#   1. PYNO.fam              — 표현형 (FID IID T1 T2 T3 T4), 6컬럼
#   2. *.grm.gz + *.grm.id   — PLINK/GCTA로 생성한 GRM
#   3. cc_cov.txt (선택)     — 분류형 공변량 (FID IID 성별 도축장 등)
#   4. qc_cov.txt (선택)     — 연속형 공변량 (FID IID PC1 PC2 등)
#
# [중요]
#   모든 파일의 FID/IID가 GRM .grm.id와 동일한 형식이어야 함
#   (예: GRM이 002136205379 이면 pheno, cc, qc도 002136205379)
#
# [사용법]
#   bash run_pipeline_v2.sh
# ═══════════════════════════════════════════════════════════════════════

set -e

# ─── 설정 (export로 Python에서도 사용 가능) ───
export GRM_PREFIX="Hanwoo_sample_96_geno_QC_final_grm_gcta"
export PYNO_FAM="PYNO.fam"
export CC_FILE="cc_cov.txt"    # 분류형 공변량 (없으면 자동 건너뜀)
export QC_FILE="qc_cov.txt"    # 연속형 공변량 (없으면 자동 건너뜀)
NUM_TRAITS=4

echo "═══════════════════════════════════════════════════════════════"
echo " MTG2 한우 4형질 다변량 유전체 예측 파이프라인 v2"
echo "═══════════════════════════════════════════════════════════════"
echo ""

# ─── 파일 존재 확인 ───
echo "[0] 입력 파일 확인..."
for f in "${GRM_PREFIX}.grm.gz" "${GRM_PREFIX}.grm.id" "${PYNO_FAM}"; do
    if [ ! -f "$f" ]; then
        echo "  ✗ $f 없음. 종료합니다."
        exit 1
    fi
    echo "  ✓ $f"
done

# 공변량 파일 확인
if [ -f "$CC_FILE" ]; then
    echo "  ✓ $CC_FILE (분류형 공변량)"
else
    echo "  - $CC_FILE 없음 (건너뜀)"
fi
if [ -f "$QC_FILE" ]; then
    echo "  ✓ $QC_FILE (연속형 공변량)"
else
    echo "  - $QC_FILE 없음 (건너뜀)"
fi
echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 1: GRM 원래 ID 기준으로 -p 파일 생성 (6컬럼 .fam 형식)
# ═══════════════════════════════════════════════════════════════
echo "[1] -p 파일 생성 (GRM 원래 ID 그대로 사용)..."

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import os

grm_prefix = os.environ['GRM_PREFIX']

# GRM ID 읽기 — 이 ID를 모든 파일의 기준으로 사용
grm_id = pd.read_csv(f"{grm_prefix}.grm.id", sep=r'\s+', header=None, dtype=str)
grm_id.columns = ['FID', 'IID']
print(f"  GRM 개체 수: {len(grm_id)}")
print(f"  ID 예시: {grm_id.iloc[0]['FID']} {grm_id.iloc[0]['IID']}")

# ★ -p 파일: 6컬럼 .fam 형식 (GRM ID 그대로)
with open('total_animals_fam.id', 'w') as f:
    for _, row in grm_id.iterrows():
        f.write(f"{row['FID']} {row['IID']} 0 0 0 -9\n")

print(f"  ✓ total_animals_fam.id (6컬럼, GRM ID 그대로)")
PYTHON_SCRIPT

echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 2: 다변량 표현형 + 공변량 파일 생성 (GRM ID 순서에 맞춤)
# ═══════════════════════════════════════════════════════════════
echo "[2] 표현형 & 공변량 파일 생성 (GRM ID 기준)..."

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import os

grm_prefix = os.environ['GRM_PREFIX']
pyno_fam = os.environ['PYNO_FAM']
cc_file = os.environ['CC_FILE']
qc_file = os.environ['QC_FILE']

def normalize_id(x):
    """ID 비교용 정규화 (앞자리 0 제거)"""
    try:
        return str(int(str(x).strip()))
    except ValueError:
        return str(x).strip()

# GRM ID 읽기 (기준)
grm_id = pd.read_csv(f"{grm_prefix}.grm.id", sep=r'\s+', header=None, dtype=str)
grm_id.columns = ['FID', 'IID']
grm_id['norm'] = grm_id['IID'].apply(normalize_id)

# ─── 표현형 매칭 ───
pyno = pd.read_csv(pyno_fam, sep=r'\s+', header=None, dtype=str)
print(f"  PYNO.fam: {len(pyno)}행, {len(pyno.columns)}컬럼")

pyno['norm'] = pyno[1].apply(normalize_id)
pyno_dict = {row['norm']: row for _, row in pyno.iterrows()}

trait_cols = [2, 3, 4, 5]
matched = 0
na_count = 0

with open('multivariate.pheno', 'w') as f:
    for _, grow in grm_id.iterrows():
        fid = grow['FID']
        iid = grow['IID']
        norm = grow['norm']

        vals = []
        if norm in pyno_dict:
            prow = pyno_dict[norm]
            has_value = False
            for ci in trait_cols:
                raw = str(prow[ci]).strip() if ci < len(prow) else 'NA'
                try:
                    fval = float(raw)
                    if pd.isna(fval):
                        vals.append('NA')
                    else:
                        vals.append(raw)
                        has_value = True
                except (ValueError, TypeError):
                    vals.append('NA')
            if has_value:
                matched += 1
            else:
                na_count += 1
        else:
            vals = ['NA'] * 4
            na_count += 1

        f.write(f"{fid} {iid} {' '.join(vals)}\n")

print(f"  multivariate.pheno: 관측={matched}, 전체NA={na_count}")

# 확인
with open('multivariate.pheno', 'r') as f:
    lines = f.readlines()
print(f"  첫 2줄:")
for l in lines[:2]:
    print(f"    {l.strip()}")
print(f"  마지막 2줄:")
for l in lines[-2:]:
    print(f"    {l.strip()}")

# ─── 공변량 매칭 (GRM ID 순서 + GRM 원래 ID로 출력) ───
for cov_path, cov_name in [(cc_file, "분류형"), (qc_file, "연속형")]:
    if not os.path.exists(cov_path):
        continue

    cov = pd.read_csv(cov_path, sep=r'\s+', header=None, dtype=str)
    cov['norm'] = cov[1].apply(normalize_id)
    cov_dict = {row['norm']: row for _, row in cov.iterrows()}

    # 공변량 컬럼 범위: 0=FID, 1=IID, 2~끝=공변량, 마지막=norm(추가열)
    cov_col_start = 2
    cov_col_end = len(cov.columns) - 1  # norm 열 직전까지

    out_path = cov_path + '.matched'
    with open(out_path, 'w') as f:
        for _, grow in grm_id.iterrows():
            fid = grow['FID']
            iid = grow['IID']
            norm = grow['norm']

            if norm in cov_dict:
                crow = cov_dict[norm]
                cov_vals = [str(crow[c]) for c in range(cov_col_start, cov_col_end)]
                f.write(f"{fid} {iid} {' '.join(cov_vals)}\n")
            else:
                cov_vals = ['NA'] * (cov_col_end - cov_col_start)
                f.write(f"{fid} {iid} {' '.join(cov_vals)}\n")

    print(f"  ✓ {out_path} ({cov_name} 공변량)")

print("  ✓ 완료")
PYTHON_SCRIPT

echo ""

# 공변량 옵션 설정 (matched 파일 사용)
COV_OPTS=""
if [ -f "${CC_FILE}.matched" ]; then
    COV_OPTS="$COV_OPTS -cc ${CC_FILE}.matched"
fi
if [ -f "${QC_FILE}.matched" ]; then
    COV_OPTS="$COV_OPTS -qc ${QC_FILE}.matched"
fi

# ═══════════════════════════════════════════════════════════════
# STEP 3: 다변량 시작값 파일 생성
# ═══════════════════════════════════════════════════════════════
echo "[3] 다변량 시작값 파일 생성..."

# -cove 1 없을 때 (잔차 공분산 미모델링)
cat > sv_multivariate.txt << 'SV_EOF'
ve 0.7
ve 0.7
ve 0.7
ve 0.7
va 0.3
va 0.3
va 0.3
va 0.3
cov 0.01
cov 0.01
cov 0.01
cov 0.01
cov 0.01
cov 0.01
SV_EOF
echo "  ✓ sv_multivariate.txt"

# -cove 1 사용 시 (잔차 공분산 포함)
cat > sv_multivariate_cove.txt << 'SV_EOF'
ve 0.7
ve 0.7
ve 0.7
ve 0.7
cov 0.01
cov 0.01
cov 0.01
cov 0.01
cov 0.01
cov 0.01
va 0.3
va 0.3
va 0.3
va 0.3
cov 0.01
cov 0.01
cov 0.01
cov 0.01
cov 0.01
cov 0.01
SV_EOF
echo "  ✓ sv_multivariate_cove.txt"
echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 4: MTG2 다변량 REML 실행
# ═══════════════════════════════════════════════════════════════
echo "[4] MTG2 다변량 REML (4형질 동시)..."
ulimit -s unlimited

echo "  실행 중... (시간이 걸릴 수 있습니다)"
echo "  명령어: mtg2 -p total_animals_fam.id -d multivariate.pheno"
echo "          -zg ${GRM_PREFIX}.grm.gz -mod 4 -cove 1 ${COV_OPTS}"
echo ""

mtg2 -p total_animals_fam.id \
     -d multivariate.pheno \
     -zg ${GRM_PREFIX}.grm.gz \
     -out multivariate_result \
     -mod 4 \
     -cove 1 \
     ${COV_OPTS} \
     -sv sv_multivariate_cove.txt 2>&1 | tee mtg2_reml_log.txt

echo ""
if [ -f "multivariate_result" ]; then
    echo "  ✓ REML 완료: multivariate_result"
    cat multivariate_result
else
    echo "  ✗ REML 실패. 로그: mtg2_reml_log.txt"
    echo ""
    echo "  [대안] -cove 없이 시도:"
    echo "  mtg2 -p total_animals_fam.id -d multivariate.pheno \\"
    echo "       -zg ${GRM_PREFIX}.grm.gz -out multivariate_result \\"
    echo "       -mod 4 ${COV_OPTS} -sv sv_multivariate.txt"
    exit 1
fi
echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 5: BLUP EBV 추출
# ═══════════════════════════════════════════════════════════════
echo "[5] BLUP EBV 추출..."

cp multivariate_result sv_converged.txt

mtg2 -p total_animals_fam.id \
     -d multivariate.pheno \
     -zg ${GRM_PREFIX}.grm.gz \
     -out multivariate_blup \
     -mod 4 \
     -cove 1 \
     ${COV_OPTS} \
     -sv sv_converged.txt \
     -bv multivariate_ebv \
     -nit 0 2>&1 | tee mtg2_blup_log.txt

if [ -f "multivariate_ebv" ]; then
    echo "  ✓ multivariate_ebv 생성 ($(wc -l < multivariate_ebv)줄)"
else
    echo "  ✗ EBV 생성 실패. 로그: mtg2_blup_log.txt"
    exit 1
fi
echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 6: 결과 정리
# ═══════════════════════════════════════════════════════════════
echo "[6] 결과 정리..."

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import os

grm_prefix = os.environ['GRM_PREFIX']
N_TRAITS = 4
trait_names = {1: '도체중', 2: '등심단면적', 3: '등지방두께', 4: '근내지방도'}

# GRM ID (원래 ID)
grm_id = pd.read_csv(f"{grm_prefix}.grm.id", sep=r'\s+', header=None, dtype=str)
grm_id.columns = ['FID', 'IID']
N = len(grm_id)

# NA 개체 찾기
na_rows = []
with open('multivariate.pheno', 'r') as f:
    for i, line in enumerate(f, 1):
        parts = line.strip().split()
        traits = parts[2:]
        if all(v == 'NA' for v in traits):
            na_rows.append(i)

print(f"  예측 대상 (전체 NA): {len(na_rows)}두")

# 다변량 EBV 파싱
ebv_data = {t: [] for t in range(1, N_TRAITS+1)}

if os.path.exists('multivariate_ebv'):
    with open('multivariate_ebv', 'r') as f:
        for line in f:
            line = line.strip()
            if not line or 'trait' in line.lower() or 'ebv' in line.lower():
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    trait_num = int(parts[0])
                    ebv_val = float(parts[1])
                    if trait_num in ebv_data:
                        ebv_data[trait_num].append(ebv_val)
                except ValueError:
                    continue

    for t in range(1, N_TRAITS+1):
        print(f"  Trait {t} ({trait_names[t]}): EBV {len(ebv_data[t])}개")

# 예측 대상 결과
results = []
for row_num in na_rows:
    idx = row_num - 1
    row = {
        '순번': row_num,
        'FID': grm_id.iloc[idx]['FID'] if idx < N else 'Unknown',
        'IID': grm_id.iloc[idx]['IID'] if idx < N else 'Unknown',
    }
    for t in range(1, N_TRAITS+1):
        col = f'{trait_names[t]}_EBV'
        row[col] = ebv_data[t][idx] if t in ebv_data and idx < len(ebv_data[t]) else None
    results.append(row)

df_14 = pd.DataFrame(results)
df_14.to_csv('prediction_results_14.csv', index=False, float_format='%.6f', encoding='utf-8-sig')

print(f"\n  ═══ 예측 대상 결과 ═══")
print(df_14.to_string(index=False, float_format='%.4f'))

# 전체 결과
with open('multivariate.pheno', 'r') as f:
    pheno_lines = f.readlines()

all_rows = []
for i in range(N):
    row = {'FID': grm_id.iloc[i]['FID'], 'IID': grm_id.iloc[i]['IID']}
    parts = pheno_lines[i].strip().split()
    traits_vals = parts[2:]
    row['표현형_유무'] = '전체NA' if all(v == 'NA' for v in traits_vals) else '있음'
    for t in range(1, N_TRAITS+1):
        col = f'{trait_names[t]}_EBV'
        row[col] = ebv_data[t][i] if t in ebv_data and i < len(ebv_data[t]) else None
    all_rows.append(row)

df_all = pd.DataFrame(all_rows)
df_all.to_csv('prediction_results_all.csv', index=False, float_format='%.6f', encoding='utf-8-sig')
print(f"\n  ✓ prediction_results_14.csv ({len(df_14)}두)")
print(f"  ✓ prediction_results_all.csv ({len(df_all)}두)")

PYTHON_SCRIPT

echo ""
echo "═══════════════════════════════════════════════════════════════"
echo " 완료!"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo " [출력 파일]"
echo "   multivariate_result       — 분산/공분산 성분"
echo "   multivariate_ebv          — 4형질 EBV (전체 개체)"
echo "   prediction_results_14.csv — 예측 대상 결과"
echo "   prediction_results_all.csv— 전체 개체 EBV"
echo "   mtg2_reml_log.txt         — REML 실행 로그"
echo "   mtg2_blup_log.txt         — BLUP 실행 로그"
echo ""
