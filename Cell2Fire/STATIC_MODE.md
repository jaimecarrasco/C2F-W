# Modo estático (`--sim-static`) — C2F-W

Rama `static-realtime`, derivada de `breaching-shading-fat`.

## Qué hace

En vez de propagar el fuego desde una ignición, recorre **todas las celdas
combustibles** del paisaje y calcula, para **cada fila de `Weather.csv`** (una
hora u escenario), las **métricas potenciales** de comportamiento del fuego:
*"si el fuego llegara a esta celda bajo este viento y humedad, así ardería"*.

Hereda del motor dinámico: el mismo kernel (`--sim S|K|C|P`), `--fmc-shading`
(humedad espacial), `--fuel-adjustment` y las tablas de lookup. No reimplementa
física: llama a `calculate_s/k/fbp`, los mismos que usa la propagación.

## Uso

```bash
./Cell2Fire --sim-static --sim S \
    --input-instance-folder <instancia>/ \
    --output-folder <salida>/ \
    --cros [--fmc-shading] [--breach-factor 2.0]
```

La instancia necesita las capas del paisaje (`fuels.asc`, `elevation`, `slope`,
`saz`, bandas de copa), la tabla de combustibles (`spain_lookup_table.csv` u
homóloga) y un `Weather.csv`.

## Salidas (una por fila de Weather.csv, sufijo = datetime)

| Capa | Significado por celda |
|---|---|
| `Intensity_<t>.asc` | intensidad de línea de Byram potencial (kW/m) |
| `ROS_<t>.asc` | velocidad de propagación de cabeza (m/min) |
| `FlameLength_<t>.asc` | longitud de llama superficial (m) |
| `MaxFlameLength_<t>.asc` | máx(superficie, copa) (m) |
| `CrownFire_<t>.asc` | 1 si iniciaría fuego de copa, 0 si no |
| `CrownFractionBurn_<t>.asc` | fracción de copa quemada [0,1] |
| `BreachRisk_<t>.asc` | 1 si la celda podría llevar fuego al otro lado de una barrera vecina (solo con `--breach-factor`) |

Celdas no-combustible: `NODATA` (-9999).

## Métricas que NO se emiten

Tiempo de llegada, dirección de avance real, celdas quemadas, trayectoria: solo
existen con propagación desde una ignición. En estático no tienen sentido.

## Riesgo de breaching

Con `--breach-factor k`, una celda combustible adyacente a una barrera (río,
camino, cortafuegos) se marca en riesgo si `k · MaxFlameLength ≥ ancho_barrera`
— el mismo criterio determinista del motor dinámico. Sirve para alertas y
recomendaciones: identifica dónde el fuego podría saltar una barrera dado el
viento de la hora.

## Tiempo real / nowcasting / escenarios

El loop por filas de `Weather.csv` no distingue pasado, presente o futuro: cada
fila es un mapa estático independiente. Para operación en tiempo real, se alimenta
el `Weather.csv` con:
- horas observadas (estación/reanálisis) → mapa actual,
- horas posteriores de un pronóstico → nowcasting,
- escenarios hipotéticos → análisis what-if.

El mapa de todo el paisaje se recalcula hora a hora con el clima de cada fila.
