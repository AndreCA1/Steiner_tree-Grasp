#!/usr/bin/env bash
# Roda o baseline (Takahashi-Matsuyama) e o GRASP em todas as instancias de benchmarks/, em paralelo.
# Uso: ./rodarTudo.sh        ou, por exemplo: TEMPO=10 EXECUCOES=3 PARALELO=8 ./rodarTudo.sh
set -euo pipefail
cd "$(dirname "$0")"

export TEMPO="${TEMPO:-60}" EXECUCOES="${EXECUCOES:-10}" SEMENTE="${SEMENTE:-1}" ALFA="${ALFA:-0.1}"
PARALELO="${PARALELO:-$(nproc)}"
export PASTA="resultados/$(date +%Y%m%d_%H%M%S)_t${TEMPO}_a${ALFA}"
mkdir -p "$PASTA"
echo "Resultados em $PASTA (tempo=$TEMPO execucoes=$EXECUCOES alfa=$ALFA paralelo=$PARALELO)"
INICIO=$(date +%s)

# Uma tarefa por (instancia, metodo); cada uma grava seu proprio CSV e log
for f in benchmarks/*/*.stp benchmarks/*/*.gr; do
    echo "$f tm"; echo "$f grasp"
done | xargs -P "$PARALELO" -L 1 sh -c '
    nome=$(basename "$0"); nome=${nome%.*}
    python3 executar.py "$0" --metodo "$1" --tempo "$TEMPO" --execucoes "$EXECUCOES" \
        --semente "$SEMENTE" --alfa "$ALFA" --csv "$PASTA/${nome}_$1.csv" > "$PASTA/${nome}_$1.txt" 2>&1 \
        && echo "ok   $nome $1" || echo "ERRO $nome $1 (veja $PASTA/${nome}_$1.txt)"
'

awk 'FNR==1 && NR!=1 {next} 1' "$PASTA"/*.csv > "$PASTA/resultados.csv"
echo "Concluido em $(( $(date +%s) - INICIO )) s -> $PASTA/resultados.csv"
