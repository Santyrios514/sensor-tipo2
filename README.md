# sensor-tipo2 · DBF 2026-27 (UPB)

Optimización de la geometría y del lastre de plomo del **sensor remolcado tipo 2**
(`modelos/analisis_tipo2.ork`): nariz + cuerpo cilíndrico + transición (*boattail*) + **tubo de
cola** cilíndrico con cuatro aletas freeform montadas sobre el tubo. Sigue la metodología de
[`dbf-sensor`](https://github.com/Santyrios514/dbf-sensor) (`sensor_opt`): barrido en malla con un
CP sustituto (Barrowman interno, sin JVM) y verificación de los mejores con OpenRocket 24.12. El
paquete es autocontenido; no importa `dbf-sensor`.

Lo que no es geometría funciona igual que en `dbf-sensor`: pared fibra de vidrio 0.3 mm + PLA
1.5 mm, electrónica (ahora de **20 mm** y 150 g) detrás del tapón de lastre delantero, herraje de remolque
de 15 g, lastre de plomo macizo con tapón delantero y trasero, amarre en el CG con tolerancia
mínima de 1.5 mm, V = 30 m/s a 1495 m (ISA) y los mismos cuatro criterios del objetivo.

## Hallazgo 1: las aletas no pueden quedar dentro del diámetro del cuerpo

El `.ork` tipo 2 tiene aletas con radio de punta igual al del cuerpo ($r_{tip} = R = 35$ mm).
OpenRocket da para ese diseño $x_{CP} = -133.6$ mm (**delante de la nariz**) y
$C_{N\alpha} = 1.25$: la transición aporta $-1.63$ y las aletas solo $+0.88$.

No es un problema de este diseño en particular. Por teoría de cuerpos esbeltos, la fuerza normal
de un tramo depende solo de la "masa aparente" de la sección al principio y al final. Una sección
cuerpo + aletas con radio de tubo $r_{tc}$ y semienvergadura $s_m$ tiene
$\pi\left(s_m^2 - r_{tc}^2 + r_{tc}^4/s_m^2\right)$, así que transición y aletas juntas aportan

$$C_{N\alpha,\,cola+aletas} = \frac{2}{R^2}\left(s_m^2 - r_{tc}^2 + \frac{r_{tc}^4}{s_m^2} - R^2\right)$$

Con $s_m = R$ vale $-2k^2(1-k^2) < 0$ para cualquier $k = r_{tc}/R < 1$: **la cola resta
sustentación detrás del CG** y desestabiliza. Barrowman llega a lo mismo: con aletas de bajo
alargamiento, $C_{N\alpha,f} \lesssim 2(1-k)^2(1+k)$, y sumado a la transición
$2(k^2-1)$ da $-2k(1-k^2) \le 0$.

`scripts/03_diagnostico_aletas.py` lo cuantifica sobre el cuerpo del `.ork` (D = 70 mm):

| $r_{tip}/R$ | D aparente [mm] | $C_{N\alpha}$ transición + aletas | $x_{CP}$ [mm] | SM máx. con plomo [cal] | masa factible [kg] |
|---|---|---|---|---|---|
| 1.00 | 70 | −0.75 | −135 | −2.8 | — |
| 1.20 | 84 | −0.16 | 26 | −0.5 | — |
| 1.45 | 101.5 | +0.70 | 135 | 1.03 | 2.1 |
| 2.00 | 140 | +2.91 | 240 | 2.49 | 6.7 |

(M = 0.09, el de 30 m/s). Con ese cuerpo hace falta $r_{tip}/R \ge 1.45$ (D aparente ≥ 101.5 mm)
solo para llegar a SM = 1, y entonces la tolerancia del amarre limita el plomo a 2.1 kg.
Ningún cuerpo de la malla es factible con $r_{tip}/R < 1.3$.

Lo que sí funciona es que las aletas no se salgan del **D acostado**: guardado con las 4 aletas a
45°, la sección ocupa un cuadrado de lado $\max(D, \sqrt2\,r_{tip})$, así que con
$r_{tip} \le D/\sqrt2$ ($r_{tip}/R \le \sqrt2 = 1.414$) el sensor ocupa en la bahía lo mismo
que el cuerpo solo. Esa es la restricción por defecto (`restricciones.D_acostado_max_rel_D: 1.0`):
el D aparente (punta a punta) queda en $\sqrt2\,D$, pero el D acostado es $D$.

## Hallazgo 2: con las aletas dentro del D acostado manda la tolerancia del amarre

En toda la frontera la restricción activa es la tolerancia del amarre, y en los diseños esbeltos
(D ≤ 80 mm) el SM queda además en su **máximo** (2.0), no con holgura. Con el amarre en el CG,

$$\text{tol}_{amarre}=\frac{\alpha_{max}\,q\,S_{ref}\,C_{N\alpha}\,SM\,D}{m\,g}\ \ge\ \text{tol}_{min}
\quad\Longrightarrow\quad
m \le \frac{\alpha_{max}\,q\,S_{ref}\,C_{N\alpha}\,SM_{max}\,D}{g\,\text{tol}_{min}} \propto D^3$$

porque $S_{ref} \propto D^2$ y $C_{N\alpha} \approx 3.0$ en toda la frontera. A D = 70 mm da
4.72 kg; el barrido da 4.74 kg. El SM "de sobra" respecto de SM = 1 es justamente el que sostiene
la tolerancia de 1.5 mm: no se puede cambiar por esbeltez. Las palancas para un sensor más
delgado con la misma masa son `SM_max_cal`, `alpha_trim_max_deg`, `tol_amarre_min_mm` y la
velocidad de remolque ($q \propto V^2$).

## Resultados de referencia (configuración por defecto, 2026-10-02)

Frontera masa–diámetro (`data_opt/frontera_masa_D.csv`): para cada D, el candidato factible de
mayor masa, siempre con D acostado = D, electrónica de 20 mm, SM ∈ [1, 2] y tolerancia de amarre
≥ 1.5 mm.

| D = D acostado [mm] | 55 | 60 | 65 | 70 | 75 | 80 | 85 | 90 |
|---|---|---|---|---|---|---|---|---|
| masa total máx. [kg] | 2.23 | 2.94 | 3.80 | 4.74 | 5.76 | 6.71 | 7.59 | 8.51 |
| SM [cal] | 2.00 | 2.00 | 2.00 | 1.99 | 1.99 | 1.97 | 1.83 | 1.74 |

El objetivo maximiza la masa hasta `masa.m_max_g` (10 kg) y **después** minimiza el D acostado.
Como 10 kg no se alcanzan, gana el D más grande de la malla. Para un sensor más esbelto, fijen
`m_max_g` en la masa que necesitan: el optimizador devuelve el D más chico que la alcanza
(`config/optimizacion_6200g.yaml`, el presupuesto con MTOW 16 kg: 6.2 kg → D = 80 mm;
4.7 kg → D = 70 mm).

Dos diseños verificados con OpenRocket 24.12 (con el plomo, la electrónica y el herraje como
*Mass components*):

| | `modelos/ganador_D90.ork` (máxima masa) | `modelos/esbelto_D70.ork` |
|---|---|---|
| Candidato | `L400_D90_n0.3_conica_Lt1.125_dtc20_tc180_m0.4_g1_s0.5_r1.41421` | `L400_D70_n0.3_conica_Lt1_dtc20_tc210_m0.55_g0.85_s0.5_r1.41421` |
| Nariz / cuerpo / transición | elipsoide 27 / 91.8 / cónica 101.2 mm (19°) | elipsoide 21 / 99 / cónica 70 mm (20°) |
| Tubo de cola | Ø 20 × 180 mm | Ø 20 × 210 mm |
| Aletas (4, Onyx 3 mm) | rectangulares 72 × 53.6 mm, r_tip = 63.6 mm | c_r 115.5, c_t 98.2, h 39.5 mm, r_tip = 49.5 mm |
| D acostado / D aparente | **90 / 127.3 mm** | **70 / 99.0 mm** |
| Masa total (plomo) | 8.51 kg (8.14 kg) | 4.74 kg (4.40 kg) |
| SM (OR) · C_Nα · C_D (OR) | 1.74 · 3.0 · 0.329 | 1.99 · 3.1 · 0.387 |

El óptimo empuja a **nariz roma, cilindro corto y un tubo de cola largo y delgado** (aletas
lejos del CG): el diámetro del tubo queda en el mínimo de 20 mm y su largo en 180–240 mm, y la
cónica gana a la ogiva y a la elipsoide del `.ork` (≈ 13 % menos masa).

## Cómo funciona

```
config/optimizacion.yaml ──► 01_optimizar_malla.py ──► data_opt/ranking.csv ──► 02_verificar_openrocket.py
                              (sustituto, sin JVM)       pareto.csv              (OpenRocket 24.12: CP, masas,
                                                                                  calibración, .ork del top 5)
                             03_diagnostico_aletas.py ──► diagnostico_aletas.csv (r_tip/R vs SM)
```

1. **Variables.** Por cuerpo: largo total $L \le 400$ mm, diámetro $D$, nariz $L_n/D$ (elipsoide),
   forma y largo de la transición $L_t/D$, diámetro $d_{tc}$ (mm) y largo $L_{tc}$ del tubo de cola. Por
   aleta: $c_r = \mu L_{tc}$, $c_t = \gamma c_r$, flecha $x_s = \sigma (c_r - c_t)$ y
   $r_{tip}/R$. La raíz termina en el extremo del tubo (como el `.ork`) y con $\sigma \le 1$ la
   punta no sobresale, así que el largo total es el del cuerpo.
2. **Perfil y cavidad.** Las funciones de forma son las de OpenRocket (incluida la transición
   `clipped`). La cavidad es la erosión morfológica del perfil por las capas de pared, y volumen y
   momento acumulados salen por Simpson:
   $$\forall(x)=\int_0^x \pi r_i^2\,dx,\qquad \Phi_1(x)=\int_0^x \pi r_i^2\,x\,dx$$
3. **CP sustituto.** Barrowman general para nariz (+2) y transición ($2(k^2-1)$) y franjas para las
   aletas, como `FinSetCalc`:
   $$C_{N\alpha,1}=\frac{2\pi s^2/A_{ref}}{1+\sqrt{1+\left(\beta s^2/(A_f\cos\Gamma)\right)^2}},\qquad
   C_{N\alpha,f}=C_{N\alpha,1}\,\frac{n}{2}\left(1+\frac{r_{tc}}{s+r_{tc}}\right)$$
   El candidato se descarta si $\sum C_{N\alpha} < \varepsilon_{CN}$ (guarda contra la división).
4. **Llenado de plomo** (idéntico a `dbf-sensor`, verificado caso por caso). Con la electrónica
   detrás del tapón delantero:
   $$x_{CG}(\ell)=\frac{M_0+\rho_b\,\Phi_1(\ell)+m_e\,\bar x_e(\ell)}{m_0+\rho_b\,\forall(\ell)+m_e},\qquad SM=\frac{x_{CP}-x_{CG}}{D}$$
   Primero el tapón delantero hasta $SM = SM_{min}$ o la geometría; si se llenó, el tapón trasero
   (transición y tubo de cola) hasta $SM = SM_{min}$; después el recorte al tope $m_{max}$ y por
   último el retroceso hasta cumplir la tolerancia del amarre
   $$\text{tol}_{amarre}=\frac{\alpha_{max}\,q\,S_{ref}\,C_{N\alpha}\,(x_{CP}-x_{CG})}{m g}\ \ge\ 1.5\text{ mm}$$
5. **Objetivo** (pesos $w = (10^6, 10^4, 10^2, 1)$, cada $f_j \in [0,1]$):
   $$J = w_1\frac{m_{max}-m}{m_{max}} + w_2\frac{D_2-D_{lo}}{D_{hi}-D_{lo}} + w_3\frac{k_{ef}-k_{lo}}{1-k_{lo}} + w_4\frac{|SM-1.5|}{0.5}$$
   con $D_2$ el D acostado $\max(D, \sqrt2\,r_{tip})$ (`objetivo.diametro: acostado`, por defecto) o
   el D aparente $2\max(R, r_{tip})$ (`aparente`), $k = d_{tc}/D$ y el $k$ efectivo de arrastre de la transición
   $k_{ef}=\sqrt{k^2+f_b(1-k^2)}$, $f_b = \mathrm{clip}\big((3 - L_t/\Delta D)/2,\,0,\,1\big)$ (Hoerner,
   como OpenRocket). SM en $[1, 2]$ es restricción dura. Con estos pesos un criterio inferior solo
   compensa ≈ 101 g de masa o ≈ 0.35 mm de D acostado (`tolerancias_implicitas.json`).
6. **Malla + refinamiento** a medio paso alrededor de los 10 mejores (≈ 1.6·10⁵ candidatos,
   ≈ 1.5 min en 4 núcleos). `frontera_masa_D.csv` guarda el mejor candidato de cada D.
7. **Verificación con OpenRocket** de los 20 mejores y una muestra estratificada de 30 por
   $(D, d_{tc})$; calibración $x_{CP}^{OR} \approx \alpha + \beta\,x_{CP}^{sust}$ y, si el residuo pasa
   de 2 mm, se recalcula la malla. El ganador es el mejor $J$ con el CP de OpenRocket.

## Validación del modelo

| Comparación | Diferencia |
|---|---|
| `.ork` base, CP total (M = 0.3) | −133.636 mm (modelo) vs −133.646 mm (OpenRocket) |
| `.ork` base, aletas | $C_{N\alpha}$ 0.88473 vs 0.88473; $x_{CP}$ 367.79 vs 367.78 mm |
| `.ork` base, masas por componente | nariz 0.7 %, el resto < 0.1 % |
| 50 candidatos verificados | CP ≤ 0.08 mm (0.07 mm tras calibrar); masa total ≤ 0.06 %; CG ≤ 0.02 mm |
| Llenado contra `dbf-sensor/sensor_opt` | 274 casos idénticos (masa, SM, CG, tolerancia, límite activo) |
| `.ork` de los ganadores reabiertos | ΔCP = 0.000 mm |

## Instalación y uso

```bash
pip install -e ".[dev]"            # numpy, scipy ≥ 1.12, pandas, pyyaml, matplotlib, pytest
pip install -e ".[openrocket]"     # script 02: orlab + JPype (JDK 17 o 21)
python -m orlab fetch 24.12        # o ORLAB_JAR=/ruta/OpenRocket-24.12.jar

python scripts/01_optimizar_malla.py                 # [--config ...] [--procesos N] [--sin-figuras]
python scripts/02_verificar_openrocket.py            # verifica, calibra y exporta los .ork
python scripts/03_diagnostico_aletas.py              # r_tip/R vs SM sobre el cuerpo del .ork
python scripts/03_diagnostico_aletas.py --cand ID    # ... o sobre un candidato del ranking
pytest                                               # los de OpenRocket se saltan sin java/orlab
```

Variantes: un YAML con `hereda: optimizacion.yaml` cambia solo lo que declara
(`config/optimizacion_6200g.yaml`: tope de 6.2 kg del presupuesto con MTOW 16 kg).

## Configuración (`config/optimizacion.yaml`)

| Sección | Qué define |
|---|---|
| `materiales`, `costos_usd_kg` | densidades (kg/m³) y precio del plomo |
| `masa` | masa TOTAL objetivo del sensor (10 kg por defecto, como la corrida de 10 kg de dbf-sensor); al alcanzarla, el optimizador busca el menor D |
| `pared` | capas de afuera hacia adentro, por estación (`nariz`, `cuerpo`, `cola`, `tubo_cola`) |
| `electronica`, `masas_puntuales` | electrónica (20 mm, 150 g) detrás del lastre y herraje de remolque |
| `lastre` | radio mínimo útil, fracción máxima de L, margen antes de la transición, tapón trasero |
| `condiciones_vuelo`, `remolque`, `envolvente` | V, altitud ISA, Mach; α de trim y tolerancia mínima del amarre; rotación de guardado |
| `geometria_fija` | nariz elipsoide, transición recortada y aletas (n, material, espesor 3 mm del .ork, flutter) |
| `malla` | listas de valores de las 11 variables |
| `restricciones` | $L_{max}$ = 400 mm, SM ∈ [1, 2], $k_{min}$, $d_{tc,min}$ = 20 mm, **D acostado ≤ D** (`D_acostado_max_rel_D`), $D_{ap,max}$ opcional, h y cuerdas mínimas, ángulo máximo de la transición, base roma |
| `objetivo` | pesos, diámetro del criterio 2 (`acostado` o `aparente`) y SM objetivo |
| `ejecucion`, `numerico`, `openrocket`, `salida` | paralelo, refinamiento, verificación, discretización, `.ork` base y carpetas |

## Salidas

`data_opt/`: `ranking.csv` (todos los candidatos, con `factible`, `motivos`, objetivos, geometría,
masas, CG, CP, SM, tolerancia de amarre y, si se verificó, los valores de OpenRocket),
`pareto.csv`, `frontera_masa_D.csv` (masa máxima por diámetro), `verificacion_or.csv`,
`calibracion.json`, `tolerancias_implicitas.json`, `diagnostico_aletas.csv` y `ork/rankNN_*.ork`
(con el plomo, la electrónica y el herraje como *Mass components*). `figs_opt/`:
`ganador_dibujo.png`, `frontera_masa_D.png`, `pareto.png`, `factibilidad.png`, `sensibilidad.png`,
`calibracion.png` y `diagnostico_aletas.png`.

## Hipótesis y limitaciones

- Barrowman/OpenRocket: subsónico, ángulos pequeños, flujo libre. No incluye la estela del avión
  ni la del cable.
- **Aletas en la estela de la transición.** Los ganadores tienen ≈ 19–20° de semiángulo
  equivalente: por encima de ~9° el flujo puede separarse y las aletas trabajar en estela, con un
  CP real más adelante que el calculado. Como las aletas quedan 100–200 mm detrás de la
  transición, el efecto es menor que con aletas pegadas a ella; si hace falta, activen
  `restricciones.angulo_cola_max_deg`.
- **El tubo de cola no se verifica estructuralmente.** Ø 20 mm con 1.8 mm de pared y 180–240 mm
  de largo, cargando las aletas: revisen rigidez a flexión y la carga de despliegue antes de
  fijarlo (si hay que engrosarlo, suban `restricciones.d_tc_min_mm` y vuelvan a correr).
- El trim es estático y lineal; la electrónica es una masa uniforme de diámetro completo.
- El plomo no se funde dentro del casco (327 °C frente a ~60–80 °C del PLA): se funde o mecaniza
  aparte y debe quedar retenido frente al despliegue y el tirón del cable.

## Estructura

```
config/     optimizacion.yaml (única fuente de parámetros), optimizacion_6200g.yaml (variante)
modelos/    analisis_tipo2.ork (geometría de referencia, OpenRocket 24.12), ganador_D90.ork y esbelto_D70.ork (verificados)
scripts/    01_optimizar_malla.py, 02_verificar_openrocket.py, 03_diagnostico_aletas.py
src/sensor_tipo2/
  config.py        YAML → SI, validación, malla de cuerpos y aletas, herencia
  perfiles.py      funciones de forma de OpenRocket; perfil nariz–cuerpo–transición–tubo
  geometria.py     erosión por capas, cavidad, masas de pared, Barrowman del cuerpo, caché por cuerpo
  aletas.py        aleta sobre el tubo de cola, Barrowman por franjas, flutter (Martin)
  lastre.py        llenado de plomo (delantero, trasero, tope de masa, tolerancia de amarre)
  sustituto.py     CP con guarda; calibración lineal contra OpenRocket
  objetivo.py      f1..f4, J, tolerancias implícitas, Pareto
  barrido.py       evaluación en paralelo, ranking, refinamiento
  verificacion.py  puente con OpenRocket (único módulo con Java), calibración, .ork
  exportar.py      CSV, JSON y figuras
tests/      pytest (geometría, aletas, lastre, objetivo, barrido e integración con OpenRocket)
```
