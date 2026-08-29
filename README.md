# Comparação de modelos de predição de siRNA: DSIR vs. OligoFormer

Este projeto compara dois modelos de predição de eficácia de siRNA — **DSIR** e **OligoFormer** —
sobre transcritos de interesse, incluindo validação contra um dataset de patente com dados
experimentais reais. É composto por dois módulos Python (um wrapper para cada modelo) e um
notebook (`compare_models.ipynb`) que os usa para gerar e comparar os resultados.

## Estrutura do repositório

```
.
├── oligoformer.py                  # wrapper para o modelo OligoFormer
├── dsir.py                         # implementação do modelo DSIR
├── compare_models.ipynb            # notebook de comparação
├── requirements.txt                # dependências Python
├── pcsk9_validation_dataset_ranked.json   # dataset de validação (patente PCSK9)
├── .gitignore
└── .env                            # variáveis de ambiente (você cria, não vai pro Git)
```

As pastas abaixo **são criadas automaticamente** na primeira execução, dentro do próprio
repositório — não é preciso criá-las manualmente nem configurar nenhum path para elas:

- `oligoformer_results/`
- `dsir_results/`
- `dsir_weights_cache/`

Recomendado adicionar essas três ao `.gitignore` — são dados gerados/cache, não código.

## Pré-requisitos

- Python 3.10 ou superior
- [conda](https://docs.conda.io/) instalado (necessário só para rodar o OligoFormer)
- Git

## Dependências Python

```
pip install -r requirements.txt
```

## Configuração

### 1. Variáveis de ambiente

Crie um arquivo `.env` na raiz do repositório com o seguinte conteúdo:

```
NCBI_EMAIL=seu-email@exemplo.com
NCBI_API_KEY=sua-api-key-aqui

OLIGOFORMER_DIR=/caminho/para/sua/instalacao/do/OligoFormer
OLIGOFORMER_CONDA_ENV=nome-do-seu-ambiente-conda
```

| Variável | Obrigatória? | Descrição |
|---|---|---|
| `NCBI_EMAIL` | Sim | E-mail de contato exigido pela API do NCBI Entrez (usada para buscar sequências) |
| `NCBI_API_KEY` | Não (recomendada) | API key do NCBI — aumenta o limite de requisições. Sem ela, o código funciona, só com um limite de requisições mais baixo |
| `OLIGOFORMER_DIR` | Sim, para rodar o OligoFormer | Diretório onde o OligoFormer foi clonado/instalado (passo 2 abaixo) |
| `OLIGOFORMER_CONDA_ENV` | Sim, para rodar o OligoFormer | Nome do ambiente conda com as dependências do OligoFormer instaladas |

Sem `NCBI_EMAIL`, o notebook falha já na primeira célula. Sem `OLIGOFORMER_DIR`/`OLIGOFORMER_CONDA_ENV`,
tudo que não depende do OligoFormer funciona normalmente — o erro só aparece na hora de efetivamente
rodar o modelo.

### 2. Instalando o OligoFormer

O OligoFormer **não faz parte deste repositório** e precisa ser instalado separadamente:

1. Clone o repositório oficial do OligoFormer.
2. Siga as instruções de setup do próprio projeto para criar o ambiente conda com as dependências
   dele.
3. Aponte `OLIGOFORMER_DIR` (passo 1) para a pasta onde ele foi clonado, e `OLIGOFORMER_CONDA_ENV`
   para o nome do ambiente conda criado.

### 3. DSIR

Não exige instalação separada. Na primeira execução, baixa os pesos do modelo do servidor oficial
do DSIR (`biodev.cea.fr`) e os guarda em cache local dentro de `dsir_weights_cache/` — é necessário
ter acesso à internet ao menos na primeira vez que cada modo (`19nt`/`21nt`) for usado.

### 4. Dataset de validação (PCSK9)

A seção de validação com a patente PCSK9, no notebook, depende do arquivo
`pcsk9_validation_dataset_ranked.json`, incluído neste repositório na raiz do projeto.

## Executando

Com o `.env` configurado e o OligoFormer instalado, abra `compare_models.ipynb` e rode as
células em ordem, de cima para baixo.

## Observações

- O OligoFormer não é determinístico: rodar o mesmo FASTA duas vezes pode gerar scores
  ligeiramente diferentes. O notebook trata isso rodando um batch de 10 execuções e reportando
  média/desvio-padrão — não é um bug do código, é uma característica do modelo.
- As células iniciais de teste de conectividade com a API do NCBI fazem chamadas reais à rede a
  cada execução do notebook; falhas transitórias (rate limit, instabilidade de rede) podem
  acontecer e normalmente basta rodar de novo.
