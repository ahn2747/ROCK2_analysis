# =============================================================
# config.R
# 프로젝트 설정 — 암종/유전자/분석 옵션은 여기서만 수정
# 새 암종 추가 / .sav 갱신 절차는 CLAUDE.md의 체크리스트 참고
# 공변량/층화/보정/제외 설정 요약은 CLAUDE.md의 "Covariates" 표 참고
# =============================================================

cfg <- list(
  cancers      = c("COAD", "READ"),        # 예: c("COAD", "READ", "LIHC")
  clinical_dir = "database",               # <CANCER>.sav 위치
  gene_dir     = "database/gene_files",    # <CANCER>_<split>_<GENE>.csv 위치
  gene_pattern = "^([A-Za-z]+)_(\\d+_\\d+)_(.+)\\.csv$",   # 암종, split, 유전자
  genes_from_gdc = c(),                    # 02a 발현 행렬에서 추가로 가져올 유전자 (예: c("CD8A"))
  primary_gene = "ROCK2",                  # 주 가설 유전자: BH 보정에서 제외, raw p 보고

  drop_columns = c("rock2g_01"),           # 중복 열 (ROCK2Group과 동일)

  # Table 1 변수 (존재하는 열만 사용)
  table1_vars = c("age", "age_g", "gender", "pathologic_stage",
                  "pathologic_t", "pathologic_n", "pathologic_m",
                  "pathologic_stage_12_34", "cea_g",
                  "lymphovascular_invasion_indicator", "vascular_invasion_indicator",
                  "residual_tumor", "histologic_diagnosis",
                  "anatomic_neoplasm_subdivision", "tumor_status"),
  table1_continuous = c("age"),            # Table 1에서 연속형으로 요약할 변수 (median [IQR], Wilcoxon)
  # 분석에서 결측(NA)으로 처리할 값 (01 recode_clinical; 03 이후 모든 분석에 적용)
  analysis_na_values = list(residual_tumor = "RX"),   # RX = 잔존 종양 평가 불가

  # ---- 공변량 / 층화 / 보정 / 제외 ---------------------------------------
  # 제외
  exclude_neoadjuvant  = FALSE,             # 신보조요법 환자 제외 (값은 schema$neoadjuvant_yes)

  # 노출(발현 그룹): 첫 번째 = 기준(reference), 두 번째 = 비교 (median 초과)
  group_levels         = c("Low", "High"),
  # factor 기준 수준 (존재하는 변수만 적용)
  reference_levels     = list(gender = "Male", pathologic_stage = "I", stage = "I",
                              pathologic_stage_12_34 = "I\u2013II"),

  # Cox 모형 (04)
  cox_covariates       = c("age", "gender", "pathologic_stage_12_34", "cea_g"),            # 다변량 보정 변수
  cox_uni_covariates   = c("age", "gender", "pathologic_stage_12_34", "cea_g"),  # 단변량 표
  covariate_scale      = list(age = list(by = 1, label = "Age (per 10 years)")),  # 연속형 단위
  stage_full           = "stage",                                 # EPV 부족 시 교체될 stage 변수
  stage_collapsed      = "pathologic_stage_12_34",                # 교체 변수 (없으면 stage로 생성)
  epv_min              = 10,               # 다변량 Cox EPV 기준 (미만 → stage 병합, 그래도 미만 → 탐색적)
  min_events           = 10,               # 유전자별 생존분석 최소 사건 수 (미만이면 건너뜀)

  # 층화 / 기관 효과
  strata_var           = "tss",            # 04 민감도 strata(), 03 TSS × 그룹 교차표
  tss_barcode_pos      = c(6, 7),          # TCGA-XX-.... 에서 기관 코드 위치
  collapse_small_levels = c(tss = 10),     # 이 수 미만 수준은 합침 (04 strata, 05 design)
  collapse_other_label = "Other",

  km_times             = c(36, 60),        # 3년/5년 생존율 (개월)
  ph_split_months      = 24,               # PH 위반 유전자: 0–24개월 / >24개월 시간 분할 Cox
  rmst_tau             = 60,               # RMST 차이 계산 시점 (개월)
  median_tie_tolerance = 1,                # median-split 재현 QC 허용 오차 (명)
  group_colors         = c(Low = "#2E6FB7", High = "#C8442F"),   # 모든 그림 공통 (group_levels 이름)

  # ---- 05 GSEA ----
  gsea_genes             = NULL,           # NULL = primary_gene, "all" = 모든 유전자, 또는 c("MS4A1", "TIMP1")
  gsea_design_covariates = c("tss"),       # 주 design: ~ <공변량> + group
  gsea_sensitivity_covariates = character(),   # 민감도 design (기본: ~ group, 공변량 없음)
  gsea_extra_collections = list(           # 기본(Hallmark, Reactome, KEGG, GO:BP) 외 추가 MSigDB 컬렉션
    C8 = list(collection = "C8")           #   C8: 세포 유형 signature (면역/B세포)
  ),
  gsea_highlight = c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_ALLOGRAFT_REJECTION",
                     "REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR",
                     "KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY"),   # 항상 enrichment plot 그릴 경로
  gsea_c8_pattern = "B_CELL|PLASMA",       # C8에서 enrichment plot을 그릴 세포 유형 (정규식, padj < 0.05)
  gsea_enrichment_top = 5,                 # enrichment plot: Hallmark 상위 N개 (padj 순) + C8 패턴 상위 N개
  gsea_size    = c(15, 500),               # fgsea minSize, maxSize
  n_cores      = 2,                        # DESeq2/fgsea 병렬 (Windows: SnowParam; 워커마다 데이터 복사 → 메모리 주의)
  seed         = 2026,

  gdc_dir       = "D:/GDCdata",            # 02a 다운로드 폴더 (Windows 경로 길이 문제 회피)
  # 02a GDC 쿼리 / 전처리 (값을 바꾸면 counts/발현 파일이 달라짐 → 05 재실행)
  gdc = list(
    project_prefix    = "TCGA-",
    data.category     = "Transcriptome Profiling",
    data.type         = "Gene Expression Quantification",
    workflow.type     = "STAR - Counts",
    count_assay       = "unstranded",      # → <C>_counts.rds (DESeq2)
    tpm_assay         = "tpm_unstrand",    # → TCGA_<C>_RNAseq_Expression.csv (log2(TPM+1))
    tumor_sample_code = "01",              # 바코드 sample type 코드 (01 = primary solid tumor)
    sample_code_pos   = c(14, 15),         # 바코드에서 sample type 위치
    patient_id_chars  = 12,                # 환자 ID = 바코드 앞 12자리
    files_per_chunk   = 20
  ),
  processed_dir = "data/processed",
  output_dir    = "output",
  manifest      = "data/processed/db_manifest.json"   # .sav 지문 (md5, n, 생존 n/사건 등)
)

# ---- 원본 열 이름 스키마 (.sav의 실제 열 이름) ------------------------------
# 모든 코드는 이 스키마로 원본 열을 찾고 표준 이름(sample_id, surv_time, status,
# pathologic_stage, age, gender, neoadjuvant, last_contact)으로 바꿔서 사용.
# 암종별로 다른 항목만 cfg$schema$<CANCER>에 적음 (modifyList로 병합).
# 선택 항목(neoadjuvant, last_contact)이 없는 암종은 NULL 대신 NA로 지정.
cfg$schema_default <- list(
  id              = "sampleID",
  time            = "Days",
  time_unit       = "days",                # "days" 또는 "months"
  status          = "Status",
  event_value     = 1,                     # 사건(사망)을 뜻하는 값 (예: 1, "Dead")
  stage           = "pathologic_Stage",
  stage_format    = "numeric",             # "numeric" (1–4) 또는 "ajcc_text" ("Stage IIIB")
  age             = "age_at_initial_pathologic_diagnosis",
  sex             = "gender",
  neoadjuvant     = "history_neoadjuvant_treatment",
  neoadjuvant_yes = "Yes",                 # 신보조요법 받음을 뜻하는 값 (대소문자 무시)
  last_contact    = "last_contact_days_to"
)
cfg$schema <- list(
  COAD = list(),
  READ = list()
  # 예: LIHC = list(stage = "ajcc_pathologic_tumor_stage", stage_format = "ajcc_text")
)

# ---- 값 라벨 / 라벨 일관성 규칙 (기본값 + 암종별 재정의) ----------------------
# 코드 → 라벨. 라벨 순서 = factor 수준 순서 (기준 수준은 reference_levels로 지정)
cfg$value_labels_default <- list(
  pathologic_stage       = c(`1` = "I", `2` = "II", `3` = "III", `4` = "IV"),
  pathologic_stage_12_34 = c(`0` = "I\u2013II", `1` = "III\u2013IV"),
  cea_g                  = c(`0` = "\u22645 ng/mL", `1` = ">5 ng/mL"),  # 5.0은 0
  age_g                  = c(Low = "\u226465", High = "\u226566")
)
# 그룹 변수(code)가 원 변수(source) > cutoff 와 일치해야 함 (불일치 시 01 중단)
# high = source > cutoff 일 때의 코드 값
cfg$label_rules_default <- list(
  cea_g                  = list(source = "cea_level_pretreatment", cutoff = 5,  high = "1"),
  age_g                  = list(source = "age",                    cutoff = 65, high = "High"),
  pathologic_stage_12_34 = list(source = "pathologic_stage",       cutoff = 2,  high = "1")
)
cfg$value_labels <- list(COAD = list(), READ = list())
cfg$label_rules  <- list(COAD = list(), READ = list())

# ---- 06 외부 검증 (GEO) ------------------------------------------------------
# 데이터셋 추가: cfg$geo_datasets에 항목 하나 추가 (코드 수정 없음)
#   fields: 표준 이름 → characteristics_ch1의 원본 key ("key: value"의 key)
#           06이 모든 원본 key를 먼저 출력하므로 그것을 보고 지정
#   필수 fields: sample_type, age, gender, stage + endpoints에 쓰는 time/event
cfg$geo_dir        <- "data/geo"               # GEOquery 다운로드 캐시 (git-ignore: data/)
cfg$geo_probe_rule <- "max_mean"               # 유전자 probe 여러 개일 때: "max_mean", "max_iqr", "max_sd"
cfg$geo_na_values  <- c("", "N/A", "NA", "ND", "na", "n/a", "NaN")   # 결측 표기

cfg$geo_datasets <- list(
  GSE39582 = list(
    platform   = "GPL570",                     # 여러 플랫폼 GSE에서 series matrix 선택 + 확인
    symbol_col = "Gene Symbol",                # GPL 주석의 유전자 기호 열 ("A /// B" 다중 매핑 허용)
    validates  = "COAD",                       # 비교할 TCGA 코호트 (output/tables/<C>/survival_summary_raw.csv)
    label      = "GSE39582 (CIT, colon)",      # 그림/표 제목
    fields = c(
      sample_type    = "dataset",              # discovery / validation / Non Tumoral
      os_time        = "os.delay (months)",
      os_event       = "os.event",
      rfs_time       = "rfs.delay",
      rfs_event      = "rfs.event",
      stage          = "tnm.stage",            # 0–4 (숫자 코드)
      age            = "age.at.diagnosis (year)",
      gender         = "Sex",
      mmr_status     = "mmr.status",           # dMMR / pMMR
      kras_mutation  = "kras.mutation",
      braf_mutation  = "braf.mutation",
      tumor_location = "tumor.location",
      adjuvant_chemo = "chemotherapy.adjuvant",
      cit_subtype    = "cit.molecularsubtype"
    ),
    exclude_sample_type = c("Non Tumoral"),    # 비종양 샘플 제외 (median split 전에)
    stage_na_values     = c("0"),              # stage 0 (4명): 모형 stage 수준에서 제외 → NA (KM/단변량에는 포함)
    # 종점: 첫 번째 = 주 종점. exclude_stage = 종점이 정의되지 않는 stage (해당 환자 제외)
    endpoints = list(
      OS  = list(time = "os_time",  event = "os_event",  time_unit = "months", event_value = 1,
                 label = "Overall survival"),
      RFS = list(time = "rfs_time", event = "rfs_event", time_unit = "months", event_value = 1,
                 label = "Relapse-free survival",
                 exclude_stage = NaN)            # stage IV: 절반이 rfs.delay = 0 (무병 상태 없음) → I–III만
    ),
    factor_levels = list(mmr_status = c("pMMR", "dMMR")),   # 첫 번째 = 기준 수준
    extra_covariate = "mmr_status",            # 다변량 + 이 변수 (TCGA에 MSI 결과가 없어서 교란 확인용)
    subgroup = list(var = "mmr_status", level = "pMMR")     # 이 하위군에서 KM/Cox 반복
  )
)
