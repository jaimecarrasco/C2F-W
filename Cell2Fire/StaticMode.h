#ifndef STATICMODE_H
#define STATICMODE_H

#include "Cells.h"
#include "ReadArgs.h"
#include <string>
#include <unordered_map>
#include <vector>

/*
 * StaticMode
 * ==========
 * Modo estatico de C2F-W (--sim-static). En lugar de propagar el fuego desde una
 * ignicion, recorre TODAS las celdas combustibles del paisaje y calcula, para
 * cada fila de Weather.csv (una hora / un escenario), las metricas POTENCIALES de
 * comportamiento del fuego: "si el fuego llegara a esta celda bajo este viento y
 * humedad, asi arderia".
 *
 * Metricas que tienen sentido por celda (intrinsecas, no requieren que el fuego
 * llegue): ROS de cabeza, intensidad de Byram, longitud de llama (superficie y
 * copa), inicio de fuego de copa, fraccion de copa quemada, llama maxima.
 *
 * Metricas que NO se emiten (solo existen con propagacion): tiempo de llegada,
 * direccion de avance real, celdas quemadas, trayectoria.
 *
 * Adicional: mapa de RIESGO DE BREACHING por celda -- una celda combustible
 * adyacente a una barrera que, dado el WD de la hora y su longitud de llama,
 * podria llevar el fuego al otro lado (k * FL >= ancho de barrera). Util para
 * alertas y recomendaciones. Reutiliza el criterio del motor dinamico.
 *
 * El loop por filas de Weather.csv admite horas observadas (tiempo real),
 * nowcasting (horas posteriores) o escenarios futuros hipoteticos: cada fila es
 * un mapa estatico independiente, sin distinguir pasado/presente/futuro.
 */

// Ejecuta el modo estatico completo. Devuelve 0 en exito.
int runStaticMode(arguments* args,
                  inputs* df,
                  int nCells,
                  fuel_coefs* coefs,
                  int rows,
                  int cols,
                  double cellSide,
                  double perimeterCells,
                  double xllcorner,
                  double yllcorner,
                  std::vector<std::vector<std::string>>& weatherData,
                  // arcos bloqueados por barrera: clave = from*nCells + to, valor = ancho (m).
                  // Breaching por arco (misma logica que el motor dinamico).
                  const std::unordered_map<long long, double>& blockedArcs,
                  const std::unordered_map<std::string, int>& fuelLookup);

#endif  // STATICMODE_H
