# ADAM Benchmark Scripts

Este diretório contém scripts para comparar a versão local do ADAM (código atual)
com a versão instalada do Bioconductor.

## O que os scripts fazem

1. Instalam e executam ADAM localmente (`--engine local`).
2. Instalam e executam ADAM do Bioconductor (`--engine bioc`).
3. Exportam:
   - tabela principal de resultados;
   - resumo de execução (tempo, número de termos);
   - `sessionInfo()`.
4. Comparam as duas saídas:
   - número de termos compartilhados;
   - sobreposição de top termos;
   - correlação de métricas numéricas.

## Uso rápido

```bash
bash scripts/benchmark/run_compare_pipeline.sh partial benchmark_results
```

Para modo completo (com bootstrap):

```bash
bash scripts/benchmark/run_compare_pipeline.sh complete benchmark_results 500
```

## Benchmark por dataset (novo)

Pipeline completo (local vs Bioconductor) para datasets externos:

```bash
bash scripts/benchmark/run_dataset_compare_pipeline.sh airway partial benchmark_results 200 3000 first
```

```bash
bash scripts/benchmark/run_dataset_compare_pipeline.sh all partial benchmark_results 200 3000 first
```

```bash
bash scripts/benchmark/run_dataset_compare_pipeline.sh sce_mock partial benchmark_results 200 3000 first
```

```bash
bash scripts/benchmark/run_dataset_compare_pipeline.sh aedes partial benchmark_results 200 3000 first
```

Parâmetros:

- `dataset`: `airway`, `all`, `sce_mock`, `aedes`
- `mode`: `partial` ou `complete`
- `bootstrap`: número de bootstraps (usado em `complete`)
- `max_genes`: limite de genes para acelerar benchmark
- `pairing`: `first` (1 par controle vs experimento) ou `all` (todos os pares)

Os scripts `run_dataset_version.R` e `compare_dataset_versions.R` também podem
ser chamados manualmente se preferir granularidade.

## Uso manual

Rodar local:

```bash
Rscript scripts/benchmark/run_adam_version.R \
  --engine local \
  --mode partial \
  --outdir benchmark_results
```

Rodar Bioconductor:

```bash
Rscript scripts/benchmark/run_adam_version.R \
  --engine bioc \
  --mode partial \
  --outdir benchmark_results
```

Comparar:

```bash
Rscript scripts/benchmark/compare_adam_versions.R \
  --mode partial \
  --outdir benchmark_results \
  --top-n 20
```

## Arquivos gerados

- `table_local_<mode>.csv`
- `table_bioc_<mode>.csv`
- `summary_local_<mode>.csv`
- `summary_bioc_<mode>.csv`
- `overview_<mode>.csv`
- `metrics_<mode>.csv`
- `session_local_<mode>.txt`
- `session_bioc_<mode>.txt`

No pipeline por dataset:

- `table_<dataset>_local_<mode>.csv`
- `table_<dataset>_bioc_<mode>.csv`
- `summary_<dataset>_local_<mode>.csv`
- `summary_<dataset>_bioc_<mode>.csv`
- `overview_<dataset>_<mode>.csv`
- `metrics_<dataset>_<mode>.csv`

## Assert automatico (PASS/FAIL)

Depois de gerar os benchmarks, rode:

```bash
Rscript scripts/benchmark/benchmark_assert.R --outdir benchmark_results
```

Arquivos:

- `benchmark_results/benchmark_assertion.csv`
- `benchmark_results/benchmark_assertion.md`

Limiares padrao:

- `shared_ratio >= 0.90`
- `top_jaccard >= 0.90`
- `spearman_mean >= 0.95` (quando aplicavel)
- `mean_abs_diff_mean <= 0.05`

## Datasets recomendados para expandir benchmark

Além de `ExpressionAedes` + `KeggPathwaysAedes`:

- RNA-seq (bulk): `airway` (SummarizedExperiment)
- Microarray: `ALL` (ExpressionSet)
- Single-cell: `TENxPBMCData::TENxPBMCData()` ou `scRNAseq::ZeiselBrainData()`
- Single-cell leve para benchmark local: `scuttle::mockSCE()` (`dataset=sce_mock`)
- Dataset interno do ADAM (offline-friendly): `ExpressionAedes` + `KeggPathwaysAedes` (`dataset=aedes`)

Observação: para datasets externos, prepare mapeamento de anotação compatível
com o ADAM (por exemplo, `DBSpecies` como arquivo/data.frame no modo `own` ou
OrgDb adequado + identificador de gene compatível).
