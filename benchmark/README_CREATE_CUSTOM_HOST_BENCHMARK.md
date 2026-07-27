# Benchmark do `Create_custom_host.sh`

Este benchmark mede a criação completa de uma referência personalizada do
MTD Explorer usando o monitor de recursos já empregado no benchmark do
instalador.

## O que é medido

- tempo total;
- CPU média e máxima da árvore de processos;
- pico de RAM e swap;
- leitura e escrita em disco;
- tráfego de rede;
- temperatura;
- tamanhos dos principais outputs;
- estado dos caches antes e depois;
- validação dos bancos Kraken2, HISAT2 e BLAST;
- validação dos recursos OrgDb, eggNOG e ssGSEA;
- linhas de erro e aviso relevantes.

## Definição da corrida

Por padrão, os outputs específicos do TaxID são apagados e reconstruídos pelo
próprio `Create_custom_host.sh`. Os downloads da espécie e os bancos
compartilhados persistentes são preservados.

O runner identifica separadamente:

- cache da referência da espécie;
- cache universal de taxonomia Kraken2;
- banco eggNOG;
- cache NCBI `gene2ensembl`.

Isso permite comparar corridas frias e quentes sem chamar ambas simplesmente
de "cold" ou "warm".

O `Create_custom_host.sh` atual usa `nproc`; portanto, o número de CPUs lógicas
detectado é registrado automaticamente.

## Exemplo para *Biomphalaria glabrata*

```bash
cd ~/MTD-Explorer

bash benchmark/run_create_custom_host_benchmark.sh \
    --machine master \
    --taxid 6526 \
    --run-number 1
```

Digite `BUILD` quando o runner mostrar os diretórios que serão reconstruídos.

## Resultados

O diretório padrão é:

```text
~/MTD_custom_host_benchmarks/
```

Cada corrida contém, entre outros:

- `custom_host_summary.txt`
- `custom_host_summary.tsv`
- `summary.tsv`
- `resource_samples.csv`
- `console.log`
- `output_validation.tsv`
- `output_inventory.tsv`
- `cache_state_before.tsv`
- `cache_state_after.tsv`
- `diagnostic_hits.tsv`
- `final_console_tail.txt`

Ao final também é criado um arquivo
`create_custom_host_benchmark_bundle_*.tar.gz`.
