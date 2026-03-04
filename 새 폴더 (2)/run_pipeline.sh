#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════
# MTG2 한우 4형질 유전체 육종가 추정 통합 파이프라인
# ═══════════════════════════════════════════════════════════════════════
# 
# [필수 입력 파일]
#   1. PYNO.fam              — 표현형 (FID IID T1 T2 T3 T4), 6컬럼
#   2. *.grm.gz + *.grm.id   — PLINK/GCTA로 생성한 GRM
#
# [사용법]
#   bash run_pipeline.sh
#
# [핵심 포인트 (삽질 기록)]
#   1. -p 파일은 반드시 6컬럼 .fam 형식 (FID IID 0 0 0 -9)
#   2. ID가 Fortran INTEGER*4 범위 초과(>2,147,483,647)하면 순번으로 재매핑
#   3. 결측치는 NA (절대 -9 사용 금지)
#   4. -mod 뒤에는 형질 수 단일 정수만 (예: -mod 1)
#   5. 단일 GRM은 -zg file.grm.gz (직접 경로)
#   6. .grm.id도 -p와 동일하게 재매핑해야 함
#   7. -sv에 숫자가 아닌 파일 경로를 줘야 EBV 생성 가능
#   8. -bv로 EBV 추출 시 -nit 0 (REML 반복 없이 BLUP만)
# ═══════════════════════════════════════════════════════════════════════

set -e  # 오류 시 즉시 중단

# ─── 설정 ───
GRM_PREFIX="Hanwoo_sample_96_geno_QC_final_grm_gcta"   # .grm.gz, .grm.id
PYNO_FAM="PYNO.fam"                                     # 표현형 파일
TRAIT_NAMES=("도체중" "등심단면적" "등지방두께" "근내지방도")
NUM_TRAITS=4

echo "═══════════════════════════════════════════════════════════════"
echo " MTG2 한우 유전체 예측 파이프라인"
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
echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 1: ID 파일 생성 (6컬럼 .fam 형식 + 순번 재매핑)
# ═══════════════════════════════════════════════════════════════
echo "[1] ID 파일 생성..."

# 원래 ID 추출 (GRM 순서 유지)
awk '{print $1, $2}' "${GRM_PREFIX}.grm.id" | tr -d '\r' > original_ids.txt
N_ANIMALS=$(wc -l < original_ids.txt)
echo "  총 개체 수: ${N_ANIMALS}"

# ★ 핵심: 순번으로 재매핑 (12자리 ID → Fortran 정수 오버플로 방지)
# ID 매핑 테이블: 원래FID 원래IID 새FID 새IID
awk '{print $1, $2, NR, NR}' original_ids.txt > id_mapping.txt

# ★ 핵심: -p 파일은 반드시 6컬럼 (FID IID 父ID 母ID 性別 表現型)
awk '{print NR, NR, 0, 0, 0, -9}' original_ids.txt > total_animals_fam.id
echo "  ✓ total_animals_fam.id (6컬럼, 순번 ID)"

# ★ 핵심: GRM의 .grm.id도 순번으로 재매핑
awk '{print NR, NR}' original_ids.txt > remapped.grm.id
ln -sf "${GRM_PREFIX}.grm.gz" remapped.grm.gz
echo "  ✓ remapped.grm.id + remapped.grm.gz (심볼릭 링크)"
echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 2: 표현형 매칭 + 개별 trait 파일 생성
# ═══════════════════════════════════════════════════════════════
echo "[2] 표현형 매칭..."

# PYNO.fam 구조: [컬럼0]=FID [컬럼1]=IID [컬럼2]=T1 [컬럼3]=T2 [컬럼4]=T3 [컬럼5]=T4
#
# ★ 핵심: GRM ID와 PYNO ID 매칭 → 순번 ID로 출력
#          결측치는 반드시 NA (MTG2가 -9를 실제값으로 인식함)
#          출력 형식: 순번FID 순번IID 형질값 (3컬럼)

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import sys

# PYNO.fam 읽기
pyno = pd.read_csv("PYNO.fam", sep=r'\s+', header=None, dtype=str)
print(f"  PYNO.fam: {len(pyno)}행, {len(pyno.columns)}컬럼")

# ID 매핑 읽기
id_map = pd.read_csv("id_mapping.txt", sep=r'\s+', header=None, dtype=str)
id_map.columns = ['orig_FID', 'orig_IID', 'new_FID', 'new_IID']

# ID 정규화 (앞자리 0 제거로 매칭률 향상)
def normalize_id(x):
    try:
        return str(int(x))
    except:
        return str(x).strip()

pyno[0] = pyno[0].apply(normalize_id)
pyno[1] = pyno[1].apply(normalize_id)
id_map['norm_orig'] = id_map['orig_IID'].apply(normalize_id)

# PYNO를 dict로 변환 (norm_id → row)
pyno_dict = {}
for _, row in pyno.iterrows():
    pyno_dict[row[1]] = row

# 형질별 파일 생성 (컬럼 2,3,4,5 = Trait 1,2,3,4)
trait_cols = [2, 3, 4, 5]
trait_names = ['도체중', '등심단면적', '등지방두께', '근내지방도']

for trait_idx, col_idx in enumerate(trait_cols):
    trait_num = trait_idx + 1
    outfile = f"trait{trait_num}_remapped.pheno"
    
    matched = 0
    na_count = 0
    
    with open(outfile, 'w') as f:
        for _, mrow in id_map.iterrows():
            norm_id = mrow['norm_orig']
            new_fid = mrow['new_FID']
            new_iid = mrow['new_IID']
            
            if norm_id in pyno_dict:
                prow = pyno_dict[norm_id]
                val = str(prow[col_idx]).strip() if col_idx < len(prow) else 'NA'
                
                # ★ 결측치 처리: 빈값, nan, -9 → 모두 NA
                try:
                    fval = float(val)
                    if pd.isna(fval) or val == '-9' or val == '-9.0':
                        f.write(f"{new_fid} {new_iid} NA\n")
                        na_count += 1
                    else:
                        f.write(f"{new_fid} {new_iid} {val}\n")
                        matched += 1
                except (ValueError, TypeError):
                    f.write(f"{new_fid} {new_iid} NA\n")
                    na_count += 1
            else:
                f.write(f"{new_fid} {new_iid} NA\n")
                na_count += 1
    
    print(f"  Trait {trait_num} ({trait_names[trait_idx]}): "
          f"관측={matched}, NA={na_count}, 합계={matched+na_count}")

print("  ✓ 표현형 매칭 완료")
PYTHON_SCRIPT

echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 3: MTG2 REML 실행 (분산성분 추정)
# ═══════════════════════════════════════════════════════════════
echo "[3] MTG2 REML 분산성분 추정..."
ulimit -s unlimited

for i in $(seq 1 $NUM_TRAITS); do
    echo ""
    echo "  ── Trait $i (${TRAIT_NAMES[$((i-1))]}) ──"
    
    # ★ 핵심 명령어:
    #   -p: 6컬럼 .fam 형식
    #   -zg: .grm.gz 직접 경로
    #   -mod 1: 형질 수 (단일 정수)
    #   -sv 0.4 0.6: 시작값 (Vg=0.4, Ve=0.6)
    
    mtg2 -p total_animals_fam.id \
         -d trait${i}_remapped.pheno \
         -zg remapped.grm.gz \
         -out trait${i}_result \
         -mod 1 \
         -sv 0.4 0.6 2>&1 | tail -5
    
    if [ -f "trait${i}_result" ]; then
        echo "  ✓ trait${i}_result 생성"
        # 분산성분 추출하여 sv 파일 생성
        grep -E "^     Ve|^     Va" trait${i}_result | awk '{print tolower($1), $2}' > sv${i}.txt
        echo "  ✓ sv${i}.txt (분산성분 저장)"
    else
        echo "  ✗ trait${i}_result 생성 실패"
    fi
done

echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 4: BLUP EBV 추출 (-nit 0)
# ═══════════════════════════════════════════════════════════════
echo "[4] BLUP EBV 추출..."

for i in $(seq 1 $NUM_TRAITS); do
    echo ""
    echo "  ── Trait $i (${TRAIT_NAMES[$((i-1))]}) ──"
    
    if [ ! -f "sv${i}.txt" ]; then
        echo "  ✗ sv${i}.txt 없음 (REML 미완료). 건너뜀."
        continue
    fi
    
    # ★ 핵심: -sv에 파일 경로, -bv로 EBV 출력, -nit 0으로 반복 없이 BLUP만
    mtg2 -p total_animals_fam.id \
         -d trait${i}_remapped.pheno \
         -zg remapped.grm.gz \
         -out trait${i}_blup \
         -mod 1 \
         -sv sv${i}.txt \
         -bv trait${i}_ebv \
         -nit 0 2>&1 | tail -3
    
    if [ -f "trait${i}_ebv" ]; then
        echo "  ✓ trait${i}_ebv 생성 ($(wc -l < trait${i}_ebv)줄)"
    else
        echo "  ✗ EBV 파일 생성 실패"
    fi
done

echo ""

# ═══════════════════════════════════════════════════════════════
# STEP 5: 결과 정리 (EBV → 원래 ID 복원, 14두 예측값 추출)
# ═══════════════════════════════════════════════════════════════
echo "[5] 결과 정리..."

python3 << 'PYTHON_SCRIPT'
import pandas as pd
import os

trait_names = {1: '도체중', 2: '등심단면적', 3: '등지방두께', 4: '근내지방도'}

# ID 매핑
id_map = pd.read_csv('id_mapping.txt', sep=r'\s+', header=None, dtype=str)
id_map.columns = ['orig_FID', 'orig_IID', 'new_FID', 'new_IID']

# NA 개체(예측 대상) 행 번호 찾기
na_rows = []
with open('trait1_remapped.pheno', 'r') as f:
    for i, line in enumerate(f, 1):
        parts = line.strip().split()
        if len(parts) >= 3 and parts[2] == 'NA':
            na_rows.append(i)

print(f"  예측 대상: {len(na_rows)}두")

# EBV 파싱 (형식: "trait번호  EBV값" 한 줄에)
ebv_data = {}
for t in range(1, 5):
    ebv_file = f'trait{t}_ebv'
    if not os.path.exists(ebv_file):
        print(f"  ✗ {ebv_file} 없음")
        continue
    
    values = []
    with open(ebv_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or 'trait' in line.lower() or 'ebv' in line.lower():
                continue
            parts = line.split()
            try:
                if len(parts) == 2:
                    values.append(float(parts[1]))
                elif len(parts) == 1:
                    values.append(float(parts[0]))
            except ValueError:
                continue
    
    ebv_data[t] = values
    print(f"  Trait {t} ({trait_names[t]}): EBV {len(values)}개 추출")

# 14두 예측 결과
results = []
for row_num in na_rows:
    idx = row_num - 1
    row = {
        '순번': row_num,
        '원래_ID': id_map.iloc[idx]['orig_IID'] if idx < len(id_map) else 'Unknown',
    }
    for t in range(1, 5):
        col = f'{trait_names[t]}_EBV'
        row[col] = ebv_data[t][idx] if t in ebv_data and idx < len(ebv_data[t]) else None
    results.append(row)

df_14 = pd.DataFrame(results)
df_14.to_csv('prediction_results_14.csv', index=False, float_format='%.6f', encoding='utf-8-sig')

print(f"\n  ═══ 14두 예측 결과 ═══")
print(df_14.to_string(index=False, float_format='%.4f'))

# 전체 결과
with open('trait1_remapped.pheno', 'r') as f:
    pheno_lines = f.readlines()

all_rows = []
for i in range(len(id_map)):
    row = {'FID': id_map.iloc[i]['orig_FID'], 'IID': id_map.iloc[i]['orig_IID']}
    parts = pheno_lines[i].strip().split()
    row['표현형_유무'] = 'NA' if parts[2] == 'NA' else '있음'
    for t in range(1, 5):
        col = f'{trait_names[t]}_EBV'
        row[col] = ebv_data[t][i] if t in ebv_data and i < len(ebv_data[t]) else None
    all_rows.append(row)

df_all = pd.DataFrame(all_rows)
df_all.to_csv('prediction_results_all.csv', index=False, float_format='%.6f', encoding='utf-8-sig')
print(f"\n  ✓ prediction_results_14.csv (14두)")
print(f"  ✓ prediction_results_all.csv ({len(df_all)}두)")

PYTHON_SCRIPT

echo ""
echo "═══════════════════════════════════════════════════════════════"
echo " 완료!"
echo "═══════════════════════════════════════════════════════════════"
echo ""
echo " [출력 파일]"
echo "   trait{1-4}_result      — 분산성분 (Ve, Va, h2)"
echo "   trait{1-4}_ebv         — 전체 개체 EBV"
echo "   prediction_results_14.csv  — 14두 예측값"
echo "   prediction_results_all.csv — 전체 6182두 EBV"
echo "   id_mapping.txt         — ID 매핑 (원래↔순번)"
echo ""
