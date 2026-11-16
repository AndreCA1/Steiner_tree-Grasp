"""Executa um metodo em uma instancia e grava uma linha por execucao no CSV.

Exemplos:
    python executar.py benchmarks/B/b01.stp --metodo tm
    python executar.py benchmarks/C/c03.stp --metodo grasp --execucoes 10 --tempo 60 --csv saida.csv
    python executar.py benchmarks/PACE/instance001.gr --metodo exato
"""
import argparse
import csv
import time
from pathlib import Path

from steiner import OTIMOS, dreyfusWagner, grasp, lerInstancia, mstPodada, takahashiMatsuyama, validar

p = argparse.ArgumentParser(description="Arvore de Steiner em grafos")
p.add_argument("instancia")
p.add_argument("--metodo", choices=["tm", "grasp", "exato"], default="grasp")
p.add_argument("--tempo", type=float, default=60, help="tempo limite por execucao do GRASP (s)")
p.add_argument("--execucoes", type=int, default=10, help="execucoes do GRASP; a execucao i usa semente + i")
p.add_argument("--semente", type=int, default=1)
p.add_argument("--alfa", type=float, default=0.1, help="tamanho da LRC do GRASP (0 = guloso, 1 = aleatorio)")
p.add_argument("--csv")
args = p.parse_args()

nome = Path(args.instancia).stem
adj, terminais = lerInstancia(args.instancia)
n, m, k = len(adj) - 1, sum(map(len, adj)) // 2, len(terminais)
ref = OTIMOS.get(nome)
print(f"{nome}: n={n} m={m} k={k} referencia={ref} metodo={args.metodo}")
if args.metodo == "exato" and k > 13:
    raise SystemExit(f"Dreyfus-Wagner e O(3^k n): inviavel para k={k} (limite 13)")

linhas = []
for i in range(args.execucoes if args.metodo == "grasp" else 1):
    semente = args.semente + i if args.metodo == "grasp" else ""
    t0 = time.perf_counter()
    iteracoes, valido = 1, True
    if args.metodo == "tm":
        custo, arvore = mstPodada(adj, takahashiMatsuyama(adj, terminais, terminais[0]), terminais)
        tempo_melhor = time.perf_counter() - t0
    elif args.metodo == "grasp":
        custo, arvore, iteracoes, tempo_melhor = grasp(adj, terminais, args.tempo, semente, args.alfa, ref)
    else:
        custo = dreyfusWagner(adj, terminais)
        tempo_melhor = time.perf_counter() - t0
    tempo = time.perf_counter() - t0
    if args.metodo != "exato":
        valido = validar(adj, terminais, arvore) == custo
    gap = 100 * (custo - ref) / ref if ref else None
    print(f"  semente={semente or '-'} custo={custo} gap={'-' if gap is None else f'{gap:.3f}%'} "
          f"valido={valido} iteracoes={iteracoes} t_melhor={tempo_melhor:.2f}s t_total={tempo:.2f}s")
    linhas.append([nome, n, m, k, args.metodo, semente, args.alfa if args.metodo == "grasp" else "", ref, custo,
                   "" if gap is None else f"{gap:.4f}", valido, iteracoes, f"{tempo_melhor:.3f}", f"{tempo:.3f}"])

if args.csv:
    novo = not Path(args.csv).exists()
    with open(args.csv, "a", newline="") as f:
        w = csv.writer(f)
        if novo:
            w.writerow(["instancia", "n", "m", "k", "metodo", "semente", "alfa", "referencia", "custo", "gap",
                        "valido", "iteracoes", "tempo_melhor", "tempo_total"])
        w.writerows(linhas)
