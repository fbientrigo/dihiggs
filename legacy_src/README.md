# Encontrando zona parametrica
Lo que se hace es buscar variando los parametros para así obtener en que parametros tenemos 

para compilar se hace un link con la libreria
```bash
g++ ParamScan.cpp -o scan -I../src -L../lib -l2HDMC -lgsl -lgslcblas -lm
```

entre las funciones utilizadas para llenar la información estna los decay Widths que relacionan la vida media y ademas el Branching Ratio

y tambien obtener la masa a partir de los parametros
```c
double THDM::get_hmass(int h) {

  if ((h<1)||(h>4)) return 0.;

  double mh[5],a,l6,l7,tb,m12_2;
  get_param_phys(mh[1],mh[2],mh[3],mh[4],a,l6,l7,m12_2,tb);
  return mh[h];
}

```

versiones que corren en paralelo
```bash
g++ ParamScanGen_pll.cpp -o paralelscan -I../src -L../lib -l2HDMC -lgsl -lgslcblas -lm -fopenmp
```

# Busqueda de Parametros
Si bien existe un rango de parametros bastante general que es permitido, la busqueda para nuestro caso en particular es acotada, donde se tiene un rango de valores más reducido.

- $\beta, m_{12}$ son cantidades cercanas a 0, en esta instancia lo dejo como $[0,0.1]$ cada una