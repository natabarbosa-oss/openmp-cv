# Coeficiente de Variação com OpenMP

Cálculo do coeficiente de variação de estatura (classes de 8 cm) e peso
(classes de 4 kg) para 1000 pessoas, com dados agrupados em intervalos de
classe, paralelizado com OpenMP.

## Compilar (MSYS2 CLANG64)
```
clang -O2 -fopenmp cv_openmp.c -o cv_par -lm
clang -O2 -Wno-unknown-pragmas cv_openmp.c -o cv_seq -lm
```

## Executar
```
./cv_par.exe dados.txt [fator_replicacao] [repeticoes]
```

## Benchmark (PowerShell)
```
powershell -ExecutionPolicy Bypass -File .\bench.ps1
```

## Arquivos
- `cv_openmp.c`: código-fonte
- `dados.txt`:    - `dados.txt`: entrada com 1000 pessoas (estatura e peso, separados por espaço). Dados sintéticos gerados pelo script abaixo.
- `gerar_dados.py`: gera o `dados.txt` (semente fixa 42; requer numpy)
- `bench.ps1` / `bench.csv`: script e resultados do benchmark
- `relatorio_openmp.pdf`: relatório