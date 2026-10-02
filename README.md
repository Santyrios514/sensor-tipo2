# sensor-tipo2 · DBF 2026-27 (UPB)

Optimización de la geometría y del lastre de plomo del **sensor remolcado tipo 2**
(`modelos/analisis_tipo2.ork`): nariz + cuerpo cilíndrico + transición (*boattail*) + **tubo de
cola** cilíndrico con cuatro aletas freeform montadas sobre el tubo. Sigue la metodología de
[`dbf-sensor`](https://github.com/Santyrios514/dbf-sensor) (`sensor_opt`): barrido en malla con un
CP sustituto (Barrowman interno, sin JVM) y verificación de los mejores con OpenRocket 24.12. El
paquete es autocontenido; no importa `dbf-sensor`.

Lo que no es geometría funciona igual que en `dbf-sensor`: pared fibra de vidrio 0.3 mm + PLA
1.5 mm, electrónica de 100 mm y 150 g detrás del tapón de lastre delantero, herraje de remolque
de 15 g, lastre de plomo macizo con tapón delantero y trasero, amarre en el CG con tolerancia
mínima de 1.5 mm, V = 30 m/s a 1495 m (ISA) y los mismos cuatro criterios del objetivo.

## Hallazgo principal: las aletas no pueden quedar dentro del diámetro del cuerpo

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
En la malla completa, **ninguno de los 50 400 candidatos con $r_{tip} \le R$ es factible**. Por
eso el optimizador deja variar $r_{tip}/R$ entre 1.0 y 2.4 y minimiza el diámetro aparente
(criterio 2), para reportar cuánto hay que salirse del diámetro del cuerpo.

## Resultados de referencia (configuración por defecto, 2026-10-02)

Ganador verificado con OpenRocket 24.12 (`modelos/ganador_10kg.ork`, con el plomo, la
electrónica y el herraje como *Mass components*):

| | |
|---|---|
| Candidato | `L400_D87.5_n0.75_conica_Lt0.7_k0.35_tc70_m1_g0.7_s1_r1.9` |
| Cuerpo | L = 400 mm, D = 87.5 mm; nariz elipsoide 65.6 mm; cuerpo 203.1 mm; transición cónica 61.25 mm (θ_eq ≈ 24.9°); tubo de cola Ø 30.6 × 70 mm |
| Aletas (4, Onyx 3 mm) | c_r = 70, c_t = 49, flecha 21, h = 67.8 mm → r_tip = 83.1 mm; **D aparente 166.3 mm**, D acostado 117.6 mm; 58 g; V_flutter ≈ 411 m/s |
| Masa | 10.000 kg, de los que 9.56 kg son plomo (8.69 kg delante, 0.87 kg detrás de la electrónica) |
| Estabilidad | x_CG = 112.6 mm, x_CP(OR) = 237.1 mm, **SM = 1.42**, C_Nα = 4.76, tolerancia de amarre 1.51 mm, C_D(OR) = 0.363 |

Qué fija el tamaño de las aletas (mismo barrido, cambiando una restricción):

| Caso | D [mm] | D aparente [mm] | r_tip/R | Nota |
|---|---|---|---|---|
| 10 kg (defecto) | 87.5 | 166 | 1.9 | tolerancia de amarre ≈ 1.5 mm, activa |
| 10 kg sin tolerancia mínima de amarre | 77.5 | 139.5 | 1.8 | SM = 1.00, tolerancia 0.62 mm |
| 10 kg con transición ≤ 15° | 90 | 180 | 2.0 | transición de 112 mm |
| 6.2 kg (`optimizacion_6200g.yaml`) | 70 | 126 | 1.8 | |

Con 10 kg, la tolerancia del amarre (≥ 1.5 mm con α ≤ 5°) exige
$C_{N\alpha}(x_{CP}-x_{CG}) \ge \text{tol}_{min}\,m g/(\alpha_{max} q S_{ref})$, que crece con la
masa: a 30 m/s el peso (98 N) domina sobre la rigidez aerodinámica
($q S_{ref} C_{N\alpha} \approx 14$ N/rad), así que más masa pide aletas más grandes.

## Cómo funciona

```
config/optimizacion.yaml ──► 01_optimizar_malla.py ──► data_opt/ranking.csv ──► 02_verificar_openrocket.py
                              (sustituto, sin JVM)       pareto.csv              (OpenRocket 24.12: CP, masas,
                                                                                  calibración, .ork del top 5)
                             03_diagnostico_aletas.py ──► diagnostico_aletas.csv (r_tip/R vs SM)
```

1. **Variables.** Por cuerpo: largo total $L \le 400$ mm, diámetro $D$, nariz $L_n/D$ (elipsoide),
   forma y largo de la transición $L_t/D$, $k = d_{tc}/D$ y largo del tubo de cola $L_{tc}$. Por
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
   $$J = w_1\frac{m_{max}-m}{m_{max}} + w_2\frac{D_{ap}-D_{lo}}{D_{hi}-D_{lo}} + w_3\frac{k_{ef}-k_{lo}}{1-k_{lo}} + w_4\frac{|SM-1.5|}{0.5}$$
   con $D_{ap} = 2\max(R, r_{tip})$ y el $k$ efectivo de arrastre de la transición
   $k_{ef}=\sqrt{k^2+f_b(1-k^2)}$, $f_b = \mathrm{clip}\big((3 - L_t/\Delta D)/2,\,0,\,1\big)$ (Hoerner,
   como OpenRocket). SM en $[1, 2]$ es restricción dura. Con estos pesos un criterio inferior solo
   compensa ≈ 101 g de masa o ≈ 1.6 mm de diámetro aparente (`tolerancias_implicitas.json`).
6. **Malla + refinamiento** a medio paso alrededor de los 10 mejores (≈ 4.2·10⁵ candidatos,
   ≈ 2.5 min en 4 núcleos).
7. **Verificación con OpenRocket** de los 20 mejores y una muestra estratificada de 30 por
   $(D, k)$; calibración $x_{CP}^{OR} \approx \alpha + \beta\,x_{CP}^{sust}$ y, si el residuo pasa
   de 2 mm, se recalcula la malla. El ganador es el mejor $J$ con el CP de OpenRocket.

## Validación del modelo

| Comparación | Diferencia |
|---|---|
| `.ork` base, CP total (M = 0.3) | −133.636 mm (modelo) vs −133.646 mm (OpenRocket) |
| `.ork` base, aletas | $C_{N\alpha}$ 0.88473 vs 0.88473; $x_{CP}$ 367.79 vs 367.78 mm |
| `.ork` base, masas por componente | nariz 0.7 %, el resto < 0.1 % |
| 50 candidatos verificados | CP ≤ 0.09 mm (0.07 mm tras calibrar); masa total ≤ 0.04 %; CG ≤ 0.03 mm |
| Llenado contra `dbf-sensor/sensor_opt` | 274 casos idénticos (masa, SM, CG, tolerancia, límite activo) |
| `.ork` del ganador reabierto | ΔCP = 0.000 mm |

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
| `masa` | tope de masa TOTAL del sensor (10 kg por defecto, como la corrida de 10 kg de dbf-sensor) |
| `pared` | capas de afuera hacia adentro, por estación (`nariz`, `cuerpo`, `cola`, `tubo_cola`) |
| `electronica`, `masas_puntuales` | electrónica detrás del lastre y herraje de remolque |
| `lastre` | radio mínimo útil, fracción máxima de L, margen antes de la transición, tapón trasero |
| `condiciones_vuelo`, `remolque`, `envolvente` | V, altitud ISA, Mach; α de trim y tolerancia mínima del amarre; rotación de guardado |
| `geometria_fija` | nariz elipsoide, transición recortada y aletas (n, material, espesor 3 mm del .ork, flutter) |
| `malla` | listas de valores de las 11 variables |
| `restricciones` | $L_{max}$ = 400 mm, SM ∈ [1, 2], $k_{min}$, $d_{tc,min}$, $D_{ap,max}$ opcional, h y cuerdas mínimas, ángulo máximo de la transición, base roma |
| `objetivo`, `ejecucion`, `numerico`, `openrocket`, `salida` | pesos, paralelo, refinamiento, verificación, discretización, `.ork` base y carpetas |

## Salidas

`data_opt/`: `ranking.csv` (todos los candidatos, con `factible`, `motivos`, objetivos, geometría,
masas, CG, CP, SM, tolerancia de amarre y, si se verificó, los valores de OpenRocket),
`pareto.csv`, `verificacion_or.csv`, `calibracion.json`, `tolerancias_implicitas.json`,
`diagnostico_aletas.csv` y `ork/rankNN_*.ork` (con el plomo, la electrónica y el herraje como
*Mass components*). `figs_opt/`: `ganador_dibujo.png`, `pareto.png`, `factibilidad.png`,
`sensibilidad.png`, `calibracion.png` y `diagnostico_aletas.png`.

## Hipótesis y limitaciones

- Barrowman/OpenRocket: subsónico, ángulos pequeños, flujo libre. No incluye la estela del avión
  ni la del cable.
- **Aletas en la estela de la transición.** Con transiciones empinadas (el ganador por defecto
  tiene ≈ 25° de semiángulo equivalente) el flujo se separa y las aletas del tubo de cola trabajan
  en estela: el CP real queda más adelante que el calculado. Revisen `theta_eq_deg` y, si hace
  falta, activen `restricciones.angulo_cola_max_deg` (p. ej. 15°).
- El trim es estático y lineal; la electrónica es una masa uniforme de diámetro completo.
- El plomo no se funde dentro del casco (327 °C frente a ~60–80 °C del PLA): se funde o mecaniza
  aparte y debe quedar retenido frente al despliegue y el tirón del cable.

## Estructura

```
config/     optimizacion.yaml (única fuente de parámetros), optimizacion_6200g.yaml (variante)
modelos/    analisis_tipo2.ork (geometría de referencia, OpenRocket 24.12) y ganador_10kg.ork (ganador verificado)
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
