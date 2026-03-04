#!/usr/bin/env python3
"""
07. MTG2 데이터 정밀 매칭 및 정규화

[수정사항 - 검토보고서 기반]
1. ❌→✅ 출력 파일명: trait2~5 → trait1~4 (WORKFLOW와 일치)
2. ❌→✅ 결측치: -9.0 → 'NA' (MTG2 매뉴얼 준수)
3. ✅ 컬럼 인덱스: PYNO.fam [0]=FID [1]=IID [2]=T1 [3]=T2 [4]=T3 [5]=T4
4. ✅ ID 정규화 함수 유지 (앞의 0 제거)
"""

import pandas as pd
import numpy as np
import os

def normalize_id(x):
    """ID 앞의 0이나 공백 차이를 완전히 제거하는 함수"""
    try:
        return str(int(float(x)))
    except:
        return str(x).strip()

def fix_all_traits():
    print("=" * 60)
    print("MTG2 데이터 정밀 매칭 및 정규화 (수정 버전)")
    print("=" * 60)

    # 1. 기준 ID 로드 (GRM .grm.id에서 생성된 total_animals.id)
    id_file = 'total_animals.id'
    if not os.path.exists(id_file):
        print(f"❌ 오류: {id_file} 파일을 찾을 수 없습니다.")
        print("  먼저 GRM .grm.id에서 ID 파일을 생성하세요:")
        print("  awk '{print $1, $2}' *.grm.id | tr -d '\\r' > total_animals.id")
        return
        
    print(f"기준 ID 파일 읽는 중: {id_file}")
    ids = pd.read_csv(id_file, sep=r'\s+', header=None, dtype=str)
    ids = ids.iloc[:, [0, 1]]
    ids.columns = ['FID', 'IID']
    ids['match_key'] = ids['IID'].apply(normalize_id)
    print(f"총 분석 대상 개체 수: {len(ids)}두")

    # 2. 원본 표현형 데이터 로드
    pheno_file = 'PYNO.fam'
    if not os.path.exists(pheno_file):
        print(f"❌ 오류: {pheno_file} 파일을 찾을 수 없습니다.")
        return

    raw_pheno = pd.read_csv(pheno_file, sep=r'\s+', header=None, dtype=str)
    print(f"표현형 파일: {pheno_file} ({raw_pheno.shape[0]}행 x {raw_pheno.shape[1]}열)")
    
    # 컬럼 구조 확인
    # PYNO.fam: [0]=FID [1]=IID [2]=T1 [3]=T2 [4]=T3 [5]=T4
    if raw_pheno.shape[1] < 6:
        print(f"⚠️ 경고: 컬럼 수가 {raw_pheno.shape[1]}개입니다. 6개(FID IID T1 T2 T3 T4) 예상.")
    
    # 3. 형질별 처리
    # ★ 수정: trait_num 1~4, col_idx 2~5
    trait_col_map = {1: 2, 2: 3, 3: 4, 4: 5}
    
    for trait_num, col_idx in trait_col_map.items():
        print(f"\n[Trait {trait_num} 처리 시작] (PYNO.fam column {col_idx})")
        
        if col_idx >= raw_pheno.shape[1]:
            print(f"  ⚠️ 컬럼 {col_idx}가 없습니다 (전체 {raw_pheno.shape[1]}컬럼)")
            continue
        
        # ID(컬럼1=IID)와 형질값(col_idx) 추출
        trait_data = raw_pheno[[1, col_idx]].copy()
        trait_data.columns = ['raw_id', 'value']
        trait_data['match_key'] = trait_data['raw_id'].apply(normalize_id)
        
        # 수치 변환
        trait_data['value'] = pd.to_numeric(trait_data['value'], errors='coerce')
        
        # -9, NA 등은 결측으로 처리
        orig_str = raw_pheno[col_idx]
        is_missing_str = orig_str.isin(['-9', 'NA', 'na', '.', ''])
        trait_data.loc[is_missing_str, 'value'] = np.nan
        
        valid = trait_data['value'].dropna()
        
        if len(valid) == 0:
            print(f"  ⚠️ Trait {trait_num}에 유효한 수치 데이터가 없습니다.")
            continue
            
        mean_v, std_v = valid.mean(), valid.std()
        print(f"  - 원본 평균: {mean_v:.2f}, 표준편차: {std_v:.2f}")
        
        # 정규화
        if std_v > 0:
            trait_data['norm_v'] = (trait_data['value'] - mean_v) / std_v
            print(f"  - 정규화 완료 (mean=0, std=1)")
        else:
            print(f"  ⚠️ 표준편차=0, 정규화 불가")
            trait_data['norm_v'] = trait_data['value']

        # 4. 기준 ID 리스트와 병합
        matched = pd.merge(ids, trait_data[['match_key', 'norm_v']], on='match_key', how='left')
        
        # 매칭 결과 진단
        success_count = matched['norm_v'].notna().sum()
        print(f"  - 매칭 성공률: {success_count} / {len(ids)} ({(success_count/len(ids))*100:.1f}%)")
        
        if success_count == 0:
            print("  ❌ 매칭 실패! ID 형식을 다시 확인하세요.")
            print(f"  [ID샘플] 기준: '{ids['match_key'].iloc[0]}', 표현형: '{trait_data['match_key'].iloc[0]}'")
            continue

        # 5. 저장
        # ★ 핵심 수정: 결측치를 NA 문자열로 코딩 (MTG2 매뉴얼 준수)
        matched['out_val'] = matched['norm_v'].apply(
            lambda x: 'NA' if pd.isna(x) else f'{x:.6f}'
        )
        
        # ★ 수정: trait1~trait4로 파일명 생성 (WORKFLOW와 일치)
        out_name = f'trait{trait_num}_matched.pheno'
        
        # MTG2 -d 형식: FID IID TRAIT_VALUE
        matched[['FID', 'IID', 'out_val']].to_csv(out_name, sep=' ', index=False, header=False)
        print(f"  ✅ {out_name} 저장 완료 (결측={len(ids)-success_count}개 → NA)")

    print("\n" + "=" * 60)
    print("✅ 모든 형질 매칭 완료!")
    print("=" * 60)
    print("\n파일 형식: FID IID TRAIT_VALUE (결측=NA)")
    print("다음 단계: bash VERIFIED_COMPLETE_WORKFLOW.sh")

if __name__ == "__main__":
    fix_all_traits()
