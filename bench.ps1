function Med($exe,$t,$f,$r) {
  $env:OMP_NUM_THREADS = $t
  $v = 1..5 | ForEach-Object {
    $l = & $exe dados.txt $f $r | Where-Object { $_ -cmatch '^TEMPO ' } | Select-Object -First 1
    [double]::Parse($l.Split(' ')[1], [Globalization.CultureInfo]::InvariantCulture)
  } | Sort-Object
  $v[2]
}
$inv = [Globalization.CultureInfo]::InvariantCulture
$out = @("N,versao,threads,tempo_s")
foreach ($c in @(@(1,20000),@(1000,20),@(10000,5))) {
  $f=$c[0]; $r=$c[1]; $N=1000*$f
  $t = (Med '.\cv_seq.exe' 1 $f $r).ToString("F9",$inv)
  $out += "$N,seq,1,$t"
  foreach ($th in 1,2,4,6,8,12) {
    $t = (Med '.\cv_par.exe' $th $f $r).ToString("F9",$inv)
    $out += "$N,par,$th,$t"
  }
  Write-Host "N=$N concluido"
}
$out | Out-File bench.csv -Encoding ascii
$out
"Pronto: bench.csv gerado"