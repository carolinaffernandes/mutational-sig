# 🧬 Assinaturas Mutacionais

Pipeline para análise de assinaturas mutacionais a partir de variantes somáticas (SNVs e InDels) usando o pacote **MutationalPatterns** no R.

------------------------------------------------------------------------

## 🔧 1. Restaurar ambiente R (renv)

Antes de começar, restaure o ambiente R para garantir a reprodutibilidade:

``` bash
Rscript -e 'renv::restore()'
```

------------------------------------------------------------------------

## 📂 2. Estrutura de Pastas

A seguinte estrutura é recomendada para o projeto:

```         
├── samples/
│   ├── variants.txt
│   ├── indels.txt
│   └── groups.txt
├── results/
│   ├── plots/
│   ├── rds/
│   └── tables/
├── snvSignatures.R
├── indelSignatures.R
└── config.yaml
```

**Descrição:** - **samples/** → arquivos de entrada (.txt) e `groups.txt`\
- **results/** → diretório de saídas\
- `plots/` → gráficos gerados com *MutationalPatterns* (customizáveis com *ggplot2*)\
- `rds/` → objetos R intermediários\
- `tables/` → resumos e conferências de contagem\
- **snvSignatures.R** → script de análise de SNVs\
- **indelSignatures.R** → script de análise de InDels\
- **config.yaml** → configurações de caminhos

------------------------------------------------------------------------

## ⚙️ 3. Arquivo de Configuração (`config.yaml`)

Exemplo de conteúdo:

``` yaml
vcf_dir: "samples"
group_file: "groups.txt"
indel_file: "indels.txt"
```

> 💡 **Observação:** Todos os arquivos devem estar no mesmo diretório indicado.

------------------------------------------------------------------------

## 📊 4. Resultados

Os resultados serão organizados automaticamente na pasta `results/`:

| Diretório        | Conteúdo                                  |
|------------------|-------------------------------------------|
| `results/plots`  | Gráficos de assinaturas mutacionais       |
| `results/rds`    | Objetos R intermediários                  |
| `results/tables` | Tabelas de resumo e contagens por amostra |

------------------------------------------------------------------------

## 📚 Referência

Blokzijl, F., Janssen, R., van Boxtel, R., & Cuppen, E. (2018).\
**MutationalPatterns: comprehensive genome-wide analysis of mutational processes.**\
*Genome Medicine, 10*(1), 33.\
<https://doi.org/10.1186/s13073-018-0539-0>

------------------------------------------------------------------------

✨ *Feito para análises reprodutíveis e visualmente elegantes com MutationalPatterns.*
