# Metodología de los planos del sensor tipo 2

Cómo se generan los planos de `planos/` (y su copia en `docs/planos/`), y qué hay que fijar para
que otra sesión de Claude Code (u otra persona) obtenga **exactamente los mismos archivos, byte a
byte**.

> **Principio.** El plano no se dibuja a mano: es una función determinista de cuatro entradas
> (código, configuración, ranking y validación) más el entorno (versiones). Si las cinco cosas son
> iguales, los PNG y PDF son idénticos y `sha256sum -c` lo comprueba. No hay que "imitar" el
> plano; hay que reproducir sus entradas.

---

## 1. Receta exacta (copiar y pegar)

```bash
git clone https://github.com/Santyrios514/sensor-tipo2 && cd sensor-tipo2
git checkout claude/sensor-tipo2-v2r2        # o el commit que acompaña a docs/planos/SHA256SUMS
python3 --version                             # 3.11.x
pip install -r requirements-lock.txt          # versiones exactas (ver §2)
python -m orlab fetch 24.12                   # OpenRocket 24.12, solo para el paso 3

rm -rf data_opt figs_opt planos               # empezar de cero (las tres están en .gitignore)
python scripts/01_optimizar_malla.py          # 1. ranking (≈ 30 s)
python scripts/03_frontera_factibilidad.py    # 2. frontera (≈ 3 min; no afecta a los planos)
python scripts/02_validar_ganadores.py \
    --cand L400_D100_n1_conica_fc0_dtc15_tc100_m1_g0.6_s1_r1.2     # 3. validación (Java)
python scripts/04_planos.py \
    --cand L400_D100_n1_conica_fc0_dtc15_tc100_m1_g0.6_s1_r1.2     # 4. planos

sha256sum -c docs/planos/SHA256SUMS_datos     # ranking, validación y perfil
(cd planos && sha256sum -c ../docs/planos/SHA256SUMS)              # planos
```

Si las dos verificaciones dicen `OK`, los planos son los mismos. Si falla la de **datos**, el
problema está antes del dibujo (código, configuración o versiones numéricas). Si solo falla la de
**planos**, el problema es del entorno gráfico (versión de matplotlib) o del paso 3 (ver §2).

**El orden importa.** El paso 3 escribe `data_opt/validacion_or.csv`, y el paso 4 lo lee para la
línea "Validación" del cajetín. Sin el paso 3 (o sin Java) los planos dicen "sin validar" y sus
huellas cambian; el resto del plano es igual. El `--cand` del paso 3 valida además la alternativa
sin banderas (puesto 137). El del paso 4 le dibuja su plano. Sin `--cand` solo salen los 5 mejores.

---

## 2. Qué fija el resultado

| Entrada | Valor de referencia | Por qué importa |
|---|---|---|
| Código | rama `claude/sensor-tipo2-v2r2`; `src/sensor_tipo2/planos.py`, `scripts/04_planos.py` | toda la geometría del dibujo está en el código |
| Configuración | `config/optimizacion.yaml` del mismo commit (D de 60 a 100 mm) | define la malla, la pared, la electrónica, el herraje, el plomo y `salida.*` |
| Ranking | `data_opt/ranking.csv` del script 01 | el plano lee la fila del candidato (CSV con `%.6g`) |
| Validación | `data_opt/validacion_or.csv` del script 02 | solo la línea "Validación" del cajetín |
| Python y librerías | Python 3.11.15; numpy 2.4.6; scipy 1.17.1; pandas 3.0.6; PyYAML 6.0.1; **matplotlib 3.11.2**; orlab 0.10.0; JPype1 1.7.1 (`requirements-lock.txt`) | matplotlib decide el rasterizado y la tipografía; numpy/scipy, los números |
| OpenRocket | `OpenRocket-24.12.jar` (sha256 `4959b72f…c4f897`), OpenJDK 21 | valores de la línea de validación |
| Tipografía | DejaVu Sans, **la que trae matplotlib** (no se usa ninguna fuente del sistema) | otra fuente cambia el ancho del texto y los bytes |
| Backend | `Agg` (forzado en `planos.py` con `matplotlib.use("Agg")`) | otro backend rasteriza distinto |
| Metadatos | PNG sin `Software`; PDF sin `CreationDate`, `Producer` ni `Creator` | sin fecha incrustada, dos corridas dan los mismos bytes |
| Resolución | `salida.dpi_planos` (no está en el YAML: 300 por defecto) | cambia el PNG, no el PDF |

No se toca `matplotlibrc`: todo lo que el dibujo necesita (tamaños, colores, grosores, tipos de
línea) se pasa explícito en el código. Un `matplotlibrc` del usuario podría cambiar los parámetros
que no se pasan (por ejemplo `hatch.linewidth`). Si las huellas no coinciden, pruebe con
`MPLCONFIGDIR=$(mktemp -d)`.

---

## 3. De la fila del ranking al dibujo

```
ranking.csv ──fila──► specs_de_fila ──► (CuerpoSpec, AletaSpec)
                                          │
                     exportar.detalle ◄───┘
                       ├─ construir_cuerpo: perfil r_e(x), cavidad r_i(x) por erosión, masas de pared
                       ├─ aletas.construir: polígono de la aleta, masa, AR
                       ├─ cp_sustituto: x_CP y C_Nα (Barrowman del modelo propio)
                       └─ lastre.llenar: ℓ (tapón delantero), ℓ₂ (trasero), CG, SM, tolerancia
                                          │
                     planos.estado ◄──────┘  qué plomo se dibuja (factible o en su SM máximo)
                                          │
                     planos.plano ───────►  rankNN_<cand_id>.png y .pdf
```

1. **Lectura de la fila** (`barrido.specs_de_fila`). Solo estas columnas definen el candidato:
   `L_mm`, `D_mm`, `L_n_rel_D`, `cola_forma`, `cola_parametro`, `f_cil`, `d_tc_mm`, `L_tc_mm`
   (cuerpo) y `mu_cr`, `gamma_ct`, `sigma_flecha`, `r_tip_rel_R` (aleta). Toda la geometría se
   **recalcula** a partir de ellas: el plano no copia cotas del CSV, y por eso T11 puede contrastar
   lo dibujado con el ranking.
2. **Cuerpo.** Malla axial uniforme de paso `numerico.dx_mm` = 0.1 mm. El perfil exterior $r_e(x)$
   sale de las funciones de forma de OpenRocket (`perfiles.py`). La cavidad $r_i(x)$ sale de la
   erosión morfológica del perfil por cada capa de pared, de afuera hacia adentro (fibra 0.3 mm +
   PLA 1.5 mm). Las fronteras entre capas (`cav.fronteras[estación][k]`) son las que se pintan.
3. **Aleta.** Polígono global (mm) con $x_{LE} = L - c_r$, $c_r = \mu L_{tc}$, $c_t = \gamma c_r$,
   $x_s = \sigma(c_r - c_t)$, $r_{tip} = (r_{tip}/R)\,R$, $h = r_{tip} - r_{tc}$, en este orden de
   vértices: $(x_{LE}, r_{tc})$, $(x_{LE}+x_s, r_{tip})$, $(x_{LE}+x_s+c_t, r_{tip})$, $(L, r_{tc})$.
4. **Lastre y electrónica** (`lastre.ModeloLastre`):
   - el tapón delantero ocupa $[x_{b0},\ x_{b0}+\ell]$, con $x_{b0}$ el primer $x$ donde
     $r_i \ge$ `lastre.r_min_util_mm` (3 mm);
   - la electrónica empieza en $x_e = \max(x_{b0}+\ell+h_e,\ x_a)$, con $h_e$ = `electronica.holgura_mm`
     (2 mm), y mide $L_e$ = 20 mm;
   - el tapón trasero ocupa $[x_{r0},\ x_{r0}+\ell_2]$, con $x_{r0} = x_e + L_e + h_e$;
   - los dos tapones llenan la cavidad completa ($-r_i$ a $+r_i$) en su tramo.
5. **Estado dibujado** (`planos.estado`):
   - **factible**: $\ell$ y $\ell_2$ del llenado; CG con los dos tapones;
   - **infactible**: $\ell$ = el que minimiza $x_{CG}$ en `numerico.n_barrido_ell` (600) puntos de
     $[0, \ell_{geo}]$, sin tapón trasero (el estado de SM máximo), y el rótulo INFACTIBLE.

---

## 4. La hoja y la escala

Todo se maqueta en **milímetros de papel** sobre una figura de matplotlib del tamaño de un A3
apaisado: `figsize = (420/25.4, 297/25.4)` pulgadas. Cada bloque es un `fig.add_axes` cuya
posición en mm se divide por (420, 297).

**Extensión de cada vista** (mm reales del sensor, con $R_{max} = \max(R, r_{tip})$):

| Vista | x | y |
|---|---|---|
| Lateral | de −34 a $L + 40$ | de $-(R_{max}+37)$ a $R_{max}+24$ |
| Posterior | de $-(R_{max}+14)$ a $R_{max}+14$ | ídem |

**Escala.** Se prueba la lista normalizada 1:1, 1:1.5, 1:2, 1:2.5, 1:3, 1:4, 1:5, 1:10 y se toma
la **primera** (la mayor) que cumpla las dos condiciones:

$$\frac{(L+74) + 2(R_{max}+14)}{e} \le 420 - 2\cdot10 - 26 = 374\ \text{mm},\qquad
\frac{\max\big(2R_{max}+61,\ 2(R_{max}+14)\big)}{e} \le 297 - 2\cdot10 - 120 = 157\ \text{mm}$$

Con $L = 400$ mm y $R_{max} = 60$ mm: ancho $(474 + 148)/e$; 1:1.5 da 415 mm (no cabe) y **1:2**
da 311 mm (cabe). Las dos vistas usan la misma escala ($e$ mm reales por mm de papel) y
`set_aspect("equal")`, así que la barra de escala vale para ambas.

**Posición de los bloques** (mm de papel, origen abajo a la izquierda):

| Bloque | x₀ | y₀ | ancho × alto |
|---|---|---|---|
| Marco | 10 | 10 | 400 × 277 (línea 1.0) |
| Vista lateral | 18 | $y_v = 297 - 10 - 14 - \max(h_{lat}, 2\,lim/e)$ | $w_{lat} = (L+74)/e$ × $h_{lat} = (2R_{max}+61)/e$ |
| Vista posterior | $18 + w_{lat} + 10$ | $y_v + (h_{lat} - 2\,lim/e)/2$ (centrada en altura) | $2\,lim/e$ × $2\,lim/e$, $lim = R_{max}+14$ |
| Títulos de las vistas | 18 y $x$ de la posterior | $y_v + \max(h_{lat}, 2\,lim/e) + 3$ | 8 pt, negrita |
| Tabla de cotas | 16 | título en $y_v - 12$, filas cada 4.6 mm desde $y_v - 18$ | 3 columnas separadas 62 mm |
| Cajetín | $420 - 10 - 205 = 205$ | 10 | 205 × 92 |
| Notas | 16 | primera línea en $10 + 92 - 4 = 98$, cada 5 mm hacia abajo | 7.5 pt el título, 6.8 pt el texto |

Títulos: "VISTA LATERAL EN CORTE (aletas en verdadera magnitud)" y "VISTA POSTERIOR (desde popa)".

---

## 5. Vista lateral en corte

Unidades de los ejes: mm reales; x desde la punta de la nariz, y = radio (arriba +, abajo −).
Se dibuja en este orden de capas (`zorder`):

| zorder | Elemento | Cómo |
|---|---|---|
| 1 | Aletas en verdadera magnitud, arriba y abajo (espejo) | `Polygon` del polígono de §3.3; relleno `#eda100` α 0.55, borde `#1f1f1e` 0.6 |
| 2 | Pared por capas, arriba y abajo | `fill_between` entre fronteras $k-1$ y $k$ de cada estación; capa 1 (exterior, fibra) `#eda100`, capa 2 (PLA) `#f5b89c`; sin borde |
| 2.5 | Tapones de plomo | `fill_between(−r_i, +r_i)` en su tramo; `#2a78d6` α 0.55, rayado `////` con borde blanco |
| 2.5 | Electrónica | `fill_between(−r_i, +r_i)` en $[x_e, x_e+L_e]$; `#1baf7a` α 0.6; texto "electrónica" vertical, 5.5 pt |
| 3 | Perfil exterior | $\pm r_e(x)$, `#1f1f1e`, 0.8 |
| 3 | Perfil interior | $\pm r_i(x)$, `#5f5e58`, 0.3 |
| 3 | Uniones (fin de nariz, inicio de transición, inicio de tubo) | segmento vertical de $-r_e$ a $+r_e$, `#5f5e58` 0.4, trazo `(0, (4, 2))` |
| 3 | Eje | de $x = -6$ a $L+6$ en y = 0, `#5f5e58` 0.4, trazo `(0, (10, 2, 2, 2))` |
| 6 | Herraje (cada `masas_puntuales`) | cuadrado negro de 4 × 4 mm centrado en $(x, 0)$; nombre con guion bajo → espacio, 5.5 pt, línea guía hasta $y = -0.55R$ |
| 7 | CG | círculo blanco (ms 6, borde 1.0) con una cruz a 45° encima; texto "CG xxx.x" 3 mm debajo |
| 7 | CP | rombo `#eb6834` (ms 5); texto "CP xxx.x" 3 mm debajo |

Colores (definidos en `exportar.py`): tinta `#1f1f1e`, tinta secundaria `#5f5e58`, rejilla
`#e4e3df`, superficie `#fcfcfb`, azul `#2a78d6`, naranja `#eb6834`, aqua `#1baf7a`, amarillo
`#eda100`; rojo de INFACTIBLE `#c8322b`.

**Barra de escala:** 5 rectángulos de 10 × 2.2 mm (reales) desde x = 0 en
$y = -(R_{max}+31)$, alternando relleno negro y blanco con borde 0.4. Rótulos "0" y "50" debajo y
"mm" a la derecha (x = 52).

---

## 6. Cotas

Todas las cotas se dibujan con dos funciones (`_Cotas.h` y `_Cotas.v`) con las mismas reglas:

- **Línea de cota:** `annotate` con `arrowstyle="<|-|>"`, `#1f1f1e`, 0.5, `mutation_scale=5`, sin
  recorte en los extremos.
- **Líneas de referencia:** desde el borde del objeto hasta 1.5 mm más allá de la línea de cota,
  `#5f5e58` 0.35.
- **Texto:** 6.5 pt, centrado, 0.8 mm por encima (horizontales) o a la izquierda/derecha
  (verticales, rotado 90°), con caja blanca sin borde (pad 0.3).
- **Formato del número** (`_txt`): dos decimales sin ceros a la derecha (`15.75`, `46.13`, `90`,
  `0`). El error máximo de lectura es 0.005 mm.
- Cada cota guarda `(valor, texto)` con su nombre para el test T11.

| Cota | Tipo | Posición | Referencias desde |
|---|---|---|---|
| `L_n`, `L_c`, `L_t`, `L_tc` (encadenadas) | horizontal | $y_1 = -(R_{max}+10)$ | $-r_e$ en cada unión (L_n también desde 0). Si $L_c = 0$: texto gris "L_c = 0 (abombado)" sobre la unión nariz–transición |
| `L` | horizontal | $y_2 = -(R_{max}+20)$ | desde la fila $y_1$ |
| `D` | vertical | x = −14, de $-R$ a $R$ | desde $x = L_n$ (fin de la nariz) |
| `d_tc` | vertical | x = $L+24$, de $-r_{tc}$ a $r_{tc}$, texto a la derecha | desde $x = L$ |
| `x_s`, `c_t` (encadenadas) | horizontal | $r_{tip}+7$ | $x_s$: desde $r_{tc}$ y $r_{tip}$; $c_t$: desde $r_{tip}$. Si $x_s = 0$ no se dibuja y vale 0 |
| `c_r` | horizontal | $r_{tip}+15$ | desde $r_{tc}$ en los dos extremos |
| `h` | vertical | x = $L+10$, de $r_{tc}$ a $r_{tip}$, texto a la derecha | desde $x = L$ y desde el borde de fuga de la punta |
| `plomo_delantero` | horizontal "plomo …" | $R+7$ | desde $r_e$ en los extremos |
| `electronica` | horizontal, texto "electrónica 20 @ x = 75" | $R+15$ | desde $r_e$ en los extremos |
| `plomo_trasero` | horizontal "plomo …" | $R+7$ | desde $r_e$ en los extremos. Si no hay, vale 0 |
| `SM_D` | horizontal, texto "SM·D = 124.9 (1.25 cal)" | y = 6, de $x_{CG}$ a $x_{CP}$ | sin referencias |
| `D_ap`, `t_aleta` | texto en la vista posterior | — | — |
| `x_electronica`, `L_electronica` | dentro del texto de `electronica` | — | — |

---

## 7. Vista posterior (desde popa)

- Círculo del cuerpo de radio $R$: relleno `#e4e3df`, borde 0.8.
- Círculo del tubo de cola de radio $r_{tc}$: relleno `#fcfcfb`, borde 0.6.
- Ejes en cruz hasta $\pm 0.9\,lim$, 0.3, trazo `(0, (10, 2, 2, 2))`.
- Cada aleta $i = 0…n-1$ con ángulo $\varphi_i = \varphi_0 + 2\pi i/n$, $\varphi_0$ =
  `geometria_fija.aletas.rotacion_deg` (45°), medido desde la vertical: rectángulo de $r_{tc}$ a
  $r_{tip}$ en la dirección $u = (\sin\varphi, \cos\varphi)$ y de espesor $t$ (3 mm) centrado.
  Relleno `#eda100`, borde 0.5.
- Círculo del D aparente $2\max(R, r_{tip})$: `#eb6834` 0.8, trazo `(0, (5, 3))`, rotulado
  "D_ap = 120" 2 mm encima.
- Pie: "4 aletas a 45° · t = 3 · r_tip = 60.0".

---

## 8. Tabla de cotas, cajetín y notas

**Tabla "COTAS [mm]"**, en este orden, repartida en 3 columnas de $\lceil n/3\rceil$ filas: L total,
L_n nariz, L_c cuerpo cilíndrico, L_t transición, L_tc tubo de cola, D cuerpo, d_tc tubo de cola,
D_ap aparente, c_r cuerda de raíz, c_t cuerda de punta, x_s flecha, h envergadura, t espesor de
aleta, tapón de plomo delantero, tapón de plomo trasero, x inicio de la electrónica, largo de la
electrónica, brazo SM·D (CG → CP). Nombre a la izquierda (6.4 pt, gris) y valor alineado a la
derecha 52 mm después. Los valores son **los mismos textos** que en las vistas.

**Cajetín** (205 × 92 mm): título en negrita (7.5 pt) sobre fondo `#e4e3df`, "Sensor remolcado tipo
2 · DBF 2026-27 (UPB) · plano de optimización", y 15 filas de igual alto (6.2 pt), con el rótulo en
gris hasta x = 0.32 y el valor desde 0.33:

| Fila | Valor | Fuente |
|---|---|---|
| Candidato | `cand_id` | fila |
| Puesto en el ranking | "1 · factible" o "— " | `puesto` |
| Masa total | `.1f` g (con "(en su SM máximo)" si es infactible) | estado |
| casco / aletas | `.1f` g / `.1f` g | `Cuerpo.m_casco`, `GeomAleta.masa` |
| plomo delantero / trasero | `.1f` g / `.1f` g | $\rho_{Pb}\forall$ de cada tapón |
| electrónica / herraje | `.1f` g / `.1f` g | configuración |
| x_CG / x_CP (desde la punta) | `.1f` / `.1f` mm | estado |
| SM · C_Nα | `.3f` cal · `.3f` | estado y `CN_alpha_total` de la fila |
| Tolerancia de amarre | `.2f` mm | llenado |
| θ_eq transición · AR aleta | `.1f`° · `.2f` | cuerpo y aleta |
| V_flutter | `.0f` m/s | fila |
| Restricción activa | `restriccion_activa` (o la falla) | fila |
| Banderas | `aletas_en_estela`, `tubo_esbelto`, `flutter_margen_bajo` que valgan True, o "ninguna" | fila |
| Validación | "OpenRocket: x_CP `.2f` mm (Δ `+.2f`), SM `.3f` (Δ `+.3f`), masa Δ `+.2f` % · VALIDADO" o "sin validar" | `validacion_or.csv` |
| Escala · unidades | "1:2 en A3 · cotas en mm" | §4 |

**Notas** (texto fijo, salvo pared, densidad y aletas, que salen de la configuración):
1. cotas en mm, x desde la punta, escala válida al imprimir el PDF en A3 al 100 %; 2. pared y plomo
rayado; 3. aletas, espesor y raíz al ras del tubo; 4. CG y CP del modelo propio, validación en el
cajetín; 5. requiere CFD o ensayo antes de fabricarse.

**Infactible:** "INFACTIBLE" en 46 pt, `#c8322b` α 0.35, rotado 12° sobre el centro de la vista
lateral, y "INFACTIBLE · falla: …" en 9 pt negrita, 4 mm debajo de la vista.

---

## 9. Salidas y nombres

| Archivo | Contenido |
|---|---|
| `planos/rankNN_<cand_id>.png` | 300 dpi (4960 × 3507 px), metadatos sin `Software` |
| `planos/rankNN_<cand_id>.pdf` | vectorial, A3, sin `CreationDate`, `Producer` ni `Creator` |
| `planos/casiNN_<cand_id>.*` | solo si no hay factibles: los 3 primeros de `casi_factibles.csv` |
| `planos/cand_<cand_id>.*` | un `--cand` infactible (sin puesto) |
| `planos/comparativo_top.png` | siluetas superpuestas, 14 × 4.2 in a 200 dpi, colores `viridis` en el orden de dibujo |
| `data_opt/rankNN_perfil.csv` | $x$, $r_e$, $r_i$ cada 0.5 mm y el polígono de la aleta, para CAD/CFD |

`NN` es el puesto con dos dígitos como mínimo (`rank01`, `rank137`). Los nombres se arman con
`ruta.parent / f"{nombre}.png"`, no con `Path.with_suffix`: el `cand_id` tiene puntos (`r1.2`) que
`with_suffix` tomaría como extensión.

---

## 10. Controles

- `pytest tests/test_scripts.py`:
  - **T11:** existen el PNG y el PDF, y cada cota dibujada coincide con su columna del ranking
    (±0.05 mm; `planos.COL_RANKING` da la correspondencia);
  - **T13:** no se escribe ningún `.ork`;
  - **`test_planos_reproducibles_byte_a_byte`:** dos dibujos del mismo candidato tienen el mismo
    sha256.
- `sha256sum -c` con las huellas de `docs/planos/` (§1).
- Revisión visual: ninguna cota tapa a otra, las aletas quedan dentro del círculo de D_ap, el
  rayado del plomo no invade la pared y el cajetín no se sale del marco.

---

## 11. Para pedírselo a otra sesión de Claude Code

> En el repo `Santyrios514/sensor-tipo2`, rama `claude/sensor-tipo2-v2r2`, sigue
> `docs/METODOLOGIA_PLANOS.md` §1 al pie de la letra: instala `requirements-lock.txt`, borra
> `data_opt/`, `figs_opt/` y `planos/`, corre los scripts 01, 03, 02 y 04 en ese orden con los
> `--cand` indicados y verifica las dos listas de huellas con `sha256sum -c`. No cambies el código
> ni la configuración. Si una huella no coincide, reporta cuál y en qué paso (datos o planos), con
> las versiones de Python, numpy, scipy, pandas y matplotlib, sin "arreglar" nada.

**Qué no hacer:**
- editar los PNG o PDF a mano;
- volver a dibujar el plano con otro código;
- cambiar fuentes, dpi o backend;
- regenerar solo el script 04 sobre un `ranking.csv` de otra configuración;
- comparar contra planos que se hicieron sin el paso 3: dicen "sin validar" y tienen otras huellas.
