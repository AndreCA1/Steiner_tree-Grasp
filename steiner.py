"""Arvore de Steiner em grafos: leitura, Takahashi-Matsuyama, GRASP, exato (Dreyfus-Wagner) e validacao.

O grafo e uma lista de dicionarios: adj[u][v] = peso da aresta (u, v), vertices de 1 a n.
Uma solucao e a adjacencia da arvore: arvore[u][v] = peso, so com os vertices usados.
"""
import heapq
import random
import time
from math import inf

# Otimos da SteinLib (steinlib.zib.de). Os do PACE foram calculados com dreyfusWagner (--metodo exato).
OTIMOS = {
    "b01": 82, "b02": 83, "b03": 138, "b04": 59,
    "c01": 85, "c02": 144, "c03": 754, "c04": 1079, "c05": 1579, "c06": 55, "c07": 102, "c08": 509,
    "d01": 106, "d02": 220, "d03": 1565, "d04": 1935, "d05": 3250, "d06": 67,
    "instance001": 503, "instance002": 111,
}


def lerInstancia(caminho):
    """Le arquivos .stp (SteinLib) e .gr (PACE 2018), que usam as mesmas linhas 'Nodes', 'E u v w' e 'T v'."""
    adj, terminais = None, []
    with open(caminho) as arquivo:
        for linha in arquivo:
            p = linha.split()
            if not p:
                continue
            if p[0] == "Nodes":
                adj = [{} for _ in range(int(p[1]) + 1)]
            elif p[0] == "E":
                u, v, w = int(p[1]), int(p[2]), int(p[3])
                if w < adj[u].get(v, inf):  # arestas paralelas: fica a mais barata
                    adj[u][v] = adj[v][u] = w
            elif p[0] == "T":
                terminais.append(int(p[1]))
    return adj, terminais


def takahashiMatsuyama(adj, terminais, inicio, rng=None, alfa=0.0):
    """Liga a arvore, um a um, ao terminal mais proximo dela pelo caminho minimo. Devolve os vertices usados.

    As distancias ate a arvore sao mantidas por um Dijkstra incremental: quando um caminho entra na arvore,
    seus vertices passam a ter distancia 0 e so eles sao reinseridos na fila.
    Com rng, o terminal e sorteado na lista restrita de candidatos (LRC) do GRASP:
    os terminais com distancia <= dmin + alfa * (dmax - dmin).
    """
    dist, pred = {inicio: 0}, {inicio: None}
    arvore, faltam = {inicio}, set(terminais) - {inicio}
    fila = [(0, inicio)]
    while True:
        while fila:
            d, u = heapq.heappop(fila)
            if d > dist[u]:
                continue
            for v, w in adj[u].items():
                if d + w < dist.get(v, inf):
                    dist[v], pred[v] = d + w, u
                    heapq.heappush(fila, (d + w, v))
        if not faltam:
            return arvore
        if rng is None:
            t = min(faltam, key=lambda x: (dist[x], x))
        else:
            dmin = min(dist[x] for x in faltam)
            dmax = max(dist[x] for x in faltam)
            t = rng.choice(sorted(x for x in faltam if dist[x] <= dmin + alfa * (dmax - dmin)))
        while t not in arvore:  # segue os predecessores ate encontrar a arvore
            arvore.add(t)
            dist[t] = 0
            heapq.heappush(fila, (0, t))
            t = pred[t]
        faltam -= arvore


def mstPodada(adj, vertices, terminais):
    """MST (Prim) do subgrafo induzido por 'vertices'; depois remove folhas que nao sao terminais.

    Devolve (custo, arvore), ou None se algum terminal ficar desconectado.
    """
    inicio = terminais[0]
    arvore = {inicio: {}}
    fila = [(w, inicio, v) for v, w in adj[inicio].items() if v in vertices]
    heapq.heapify(fila)
    while fila:
        w, u, v = heapq.heappop(fila)
        if v in arvore:
            continue
        arvore[v] = {u: w}
        arvore[u][v] = w
        for x, wx in adj[v].items():
            if x in vertices and x not in arvore:
                heapq.heappush(fila, (wx, v, x))
    if any(t not in arvore for t in terminais):
        return None

    eh_terminal = set(terminais)
    folhas = [v for v in arvore if len(arvore[v]) == 1 and v not in eh_terminal]
    while folhas:
        v = folhas.pop()
        (u, _), = arvore.pop(v).items()
        del arvore[u][v]
        if len(arvore[u]) == 1 and u not in eh_terminal:
            folhas.append(u)
    return custoArvore(arvore), arvore


def custoArvore(arvore):
    return sum(w for u in arvore for w in arvore[u].values()) // 2


def buscaLocal(adj, terminais, solucao, prazo, rng=None):
    """Insercao e remocao de vertices de Steiner com primeira melhora; cada vizinho e avaliado por mstPodada.

    Candidatos: vertices de Steiner da arvore (remocao) e vertices fora dela com >= 2 vizinhos na arvore
    (insercao; com menos vizinhos o vertice seria podado). Repete enquanto houver melhora e tempo.
    """
    custo, arvore = solucao
    eh_terminal = set(terminais)
    melhorou = True
    while melhorou:
        melhorou = False
        vertices = set(arvore)
        candidatos = [v for v in vertices if v not in eh_terminal]
        candidatos += list({x for u in vertices for x in adj[u] if x not in vertices})
        if rng:
            rng.shuffle(candidatos)
        for v in candidatos:
            if time.perf_counter() > prazo:
                return custo, arvore
            if v in vertices:
                teste = vertices - {v}
            elif sum(1 for x in adj[v] if x in vertices) >= 2:
                teste = vertices | {v}
            else:
                continue
            nova = mstPodada(adj, teste, terminais)
            if nova and nova[0] < custo:
                custo, arvore = nova
                vertices = set(arvore)
                melhorou = True
    return custo, arvore


def grasp(adj, terminais, tempo, semente, alfa, alvo=None):
    """GRASP: Takahashi-Matsuyama aleatorizado (inicio sorteado + LRC) seguido de busca local.

    Para no tempo limite ou ao atingir 'alvo' (o otimo conhecido, que nao pode ser melhorado).
    Devolve (custo, arvore, iteracoes, tempo_do_melhor).
    """
    rng = random.Random(semente)
    t0 = time.perf_counter()
    prazo = t0 + tempo
    melhor, iteracoes, tempo_melhor = None, 0, 0.0
    while time.perf_counter() < prazo and not (melhor and alvo and melhor[0] <= alvo):
        vertices = takahashiMatsuyama(adj, terminais, rng.choice(terminais), rng, alfa)
        solucao = buscaLocal(adj, terminais, mstPodada(adj, vertices, terminais), prazo, rng)
        iteracoes += 1
        if melhor is None or solucao[0] < melhor[0]:
            melhor, tempo_melhor = solucao, time.perf_counter() - t0
    return melhor[0], melhor[1], iteracoes, tempo_melhor


def dreyfusWagner(adj, terminais):
    """Custo otimo por programacao dinamica, O(3^k n + 2^k m log n). Viavel so para poucos terminais (k <= ~10).

    custo[S][v] = menor arvore que conecta os terminais do conjunto S e o vertice v.
    """
    n, k = len(adj) - 1, len(terminais)

    def dijkstra(d):
        fila = [(d[v], v) for v in range(1, n + 1) if d[v] < inf]
        heapq.heapify(fila)
        while fila:
            du, u = heapq.heappop(fila)
            if du > d[u]:
                continue
            for v, w in adj[u].items():
                if du + w < d[v]:
                    d[v] = du + w
                    heapq.heappush(fila, (d[v], v))
        return d

    custo = {}
    for S in range(1, 1 << k):  # ordem crescente: os subconjuntos de S ja foram calculados
        d = [inf] * (n + 1)
        if S & (S - 1) == 0:
            d[terminais[S.bit_length() - 1]] = 0
        else:
            sub = (S - 1) & S
            while sub:
                if sub < S ^ sub:  # cada particao {sub, S - sub} uma vez so
                    a, b = custo[sub], custo[S ^ sub]
                    for v in range(1, n + 1):
                        if a[v] + b[v] < d[v]:
                            d[v] = a[v] + b[v]
                sub = (sub - 1) & S
        custo[S] = dijkstra(d)
    return min(custo[(1 << k) - 1])


def validar(adj, terminais, arvore):
    """Confere, de forma independente, se 'arvore' e uma arvore de Steiner valida. Devolve o custo ou None."""
    for u in arvore:
        for v, w in arvore[u].items():
            if adj[u].get(v) != w or arvore.get(v, {}).get(u) != w:
                return None  # aresta inexistente no grafo, peso errado ou adjacencia assimetrica
    arestas = sum(len(viz) for viz in arvore.values()) // 2
    visitados, pilha = {terminais[0]}, [terminais[0]]
    while pilha:
        for v in arvore.get(pilha.pop(), {}):
            if v not in visitados:
                visitados.add(v)
                pilha.append(v)
    conexa = len(visitados) == len(arvore)
    if not conexa or arestas != len(arvore) - 1 or any(t not in visitados for t in terminais):
        return None
    return custoArvore(arvore)
