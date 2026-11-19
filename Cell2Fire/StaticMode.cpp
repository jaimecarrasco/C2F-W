#include "StaticMode.h"
#include "ReadCSV.h"
#include "ScottAndBurganKernel.h"
#include "CanadianFBPKernel.h"
#include "KitralKernel.h"
#include "FuelModelUtils.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// Forward decls de los kernels (mismas firmas que usa Cells.cpp)
void calculate_s(inputs*, fuel_coefs*, arguments*, main_outs*, snd_outs*, fire_struc*,
                 fire_struc*, fire_struc*, bool&, weatherDF*);
void calculate_k(inputs*, inputs*, int, fuel_coefs*, arguments*, main_outs*, snd_outs*,
                 fire_struc*, fire_struc*, fire_struc*, bool&, weatherDF*);
void calculate_fbp(inputs*, fuel_coefs*, main_outs*, snd_outs*, fire_struc*, fire_struc*,
                   fire_struc*, weatherDF*);

static const char* NF = "NF";
static const char* ND = "ND";

// ---- Fisica direccional (elipse de fuego), replicada de Cells.cpp ----
// rhoTheta: radio de la elipse a un angulo 'theta' (grados) desde el eje mayor,
// con semiejes a (mayor) y b (menor). Da el ROS en esa direccion.
static double rhoTheta(double theta, double a, double b)
{
    if (a <= 0.0 || b <= 0.0 || a < b)
        return (a > 0.0 ? a : 0.0);
    const double pi = 3.141592653589793;
    double c2 = a * a - b * b;
    double e = std::sqrt(c2) / a;
    double r1 = a * (1.0 - e * e);
    double r2 = 1.0 - e * std::cos(theta * pi / 180.0);
    return (r2 != 0.0) ? r1 / r2 : a;
}

// Angulo (convencion de la grilla, 0=Este,90=Norte,...) hacia cada uno de los 8
// vecinos, indexado por (dr,dc). Coincide con angleDict de Cells.cpp.
static int neighborAngle(int dr, int dc)
{
    // dr: fila (hacia abajo positivo), dc: columna (hacia derecha positivo)
    if (dr == 0 && dc == 1) return 0;     // E
    if (dr == -1 && dc == 1) return 45;   // NE
    if (dr == -1 && dc == 0) return 90;   // N
    if (dr == -1 && dc == -1) return 135; // NO
    if (dr == 0 && dc == -1) return 180;  // O
    if (dr == 1 && dc == -1) return 225;  // SO
    if (dr == 1 && dc == 0) return 270;   // S
    if (dr == 1 && dc == 1) return 315;   // SE
    return -1;
}

// Direccion de cabeza en convencion cartesiana de la grilla, a partir del azimut
// meteorologico de avance (waz/raz). Replica el bloque de Cells.cpp (lineas ~553-567).
static double headingCartesian(double waz)
{
    double offset = waz + 270.0;
    double ca = 360.0 - (offset >= 360.0 ? (waz - 90.0) : offset);
    if (ca == 360.0) ca = 0.0;
    if (ca < 0.0) ca += 360.0;
    return ca;
}

// Escribe una capa raster ascii-grid con la cabecera estandar.
static void writeAsc(const std::string& path, const std::vector<float>& data, int rows,
                     int cols, double xll, double yll, double cell, float nodata = -9999.0f)
{
    std::ofstream f(path);
    f << "ncols " << cols << "\n";
    f << "nrows " << rows << "\n";
    f << "xllcorner " << std::fixed << std::setprecision(6) << xll << "\n";
    f << "yllcorner " << yll << "\n";
    f << "cellsize " << cell << "\n";
    f << "NODATA_value " << (int)nodata << "\n";
    f << std::setprecision(4);
    for (int r = 0; r < rows; ++r)
    {
        for (int c = 0; c < cols; ++c)
        {
            float v = data[r * cols + c];
            f << v;
            f << (c + 1 < cols ? ' ' : '\n');
        }
    }
}

// Etiqueta de hora segura para nombre de archivo (datetime -> YYYYMMDD_HHMM o indice).
static std::string hourTag(const std::vector<std::string>& row, int idx)
{
    std::string dt = (row.size() > 1 ? row[1] : "");
    std::string clean;
    for (char ch : dt)
    {
        if (std::isalnum((unsigned char)ch))
            clean += ch;
        else if (ch == ' ' || ch == '-' || ch == ':' || ch == '/')
            clean += (clean.empty() || clean.back() == '_') ? "" : "_";
    }
    if (clean.empty())
    {
        std::ostringstream os;
        os << "t" << std::setw(3) << std::setfill('0') << idx;
        clean = os.str();
    }
    return clean;
}

int
runStaticMode(arguments* args, inputs* df, int nCells, fuel_coefs* coefs, int rows,
              int cols, double cellSide, double perimeterCells, double xllcorner, double yllcorner,
              std::vector<std::vector<std::string>>& weatherData,
              const std::unordered_map<long long, double>& blockedArcs,
              const std::unordered_map<std::string, int>& /*fuelLookup*/)
{
    // Cargar Weather.csv directamente: en modo estatico no corre el loop dinamico
    // que normalmente puebla WeatherData, asi que lo leemos aqui.
    std::string wpath = args->InFolder;
    if (!wpath.empty() && wpath.back() != '/')
        wpath += "/";
    wpath += "Weather.csv";
    CSVReader reader(wpath, ",");
    weatherData = reader.getData(wpath);

    const int WPeriods = static_cast<int>(weatherData.size()) - 1;  // fila 0 = header
    if (WPeriods < 1)
    {
        std::cerr << "[sim-static] Weather.csv sin filas de datos (" << wpath << ").\n";
        return 1;
    }

    // Parsear TODAS las filas de clima una sola vez (reutiliza el parser del motor).
    std::vector<weatherDF> wdf;
    reader.parseWeatherDF(wdf, args, weatherData, WPeriods);

    std::string outDir = args->OutFolder;
    if (!outDir.empty() && outDir.back() != '/')
        outDir += "/";

    std::cout << "[sim-static] " << WPeriods << " escenarios (filas de Weather.csv), "
              << nCells << " celdas, kernel " << args->Simulator << "\n";

    // Buffers por capa (reutilizados cada hora)
    std::vector<float> rosL(nCells), intenL(nCells), flL(nCells), crownL(nCells),
        cfbL(nCells), maxflL(nCells), breachL(nCells);

    // Grafo de aristas (ROS direccional) para rutas de fuego minimas, si se pidio.
    // Formato: una fila por arista i->j con el ROS de i hacia j y el tiempo de
    // transito (distancia centro-a-centro / ROS). Peso recomendado para Dijkstra:
    // la columna 'time_min'. Se emite un CSV por escenario.
    const double diagDist = cellSide * std::sqrt(2.0);

    // Loop temporal: una capa-conjunto por fila de clima (hora observada / nowcast /
    // escenario futuro -- cada fila es un mapa estatico independiente).
    for (int wp = 0; wp < WPeriods; ++wp)
    {
        weatherDF* wptr = &wdf[wp];

        std::ofstream edgeFile;
        if (args->StaticGraph)
        {
            std::string tag = hourTag(weatherData[wp + 1], wp);
            edgeFile.open(outDir + "FireGraph_" + tag + ".csv");
            // i->j; ros direccional (m/min); distancia (m); tiempo de transito (min);
            // intensidad de Byram direccional (kW/m, escalada por ROS_dir/ROS_cabeza).
            // Pesos para NetworkX: 'time_min' (rutas mas rapidas) o 'fli_kw_m' (max intensidad).
            edgeFile << "i,j,ros_m_min,dist_m,time_min,fli_kw_m\n";
        }

        for (int cid = 1; cid <= nCells; ++cid)
        {
            inputs& d = df[cid - 1];
            // No-combustible: NODATA en todas las capas. Se filtra por el string de
            // tipo de combustible (robusto a la numeracion interna, que para FBP marca
            // NF con ftypeN=101). Captura "NF"/"ND" por prefijo y nftype no positivo.
            if (std::strncmp(d.fueltype, "NF", 2) == 0 || std::strncmp(d.fueltype, "ND", 2) == 0 ||
                d.fueltype[0] == '\0' || d.nftype <= 0)
            {
                rosL[cid - 1] = intenL[cid - 1] = flL[cid - 1] = -9999.0f;
                crownL[cid - 1] = cfbL[cid - 1] = maxflL[cid - 1] = -9999.0f;
                breachL[cid - 1] = -9999.0f;
                continue;
            }

            // El viento de la celda proviene del escenario de la hora (campo uniforme;
            // con --fmc-shading la humedad se acondiciona por celda dentro del kernel).
            d.ws = wptr->ws;
            d.waz = (int)wptr->waz;

            main_outs mainstruct = {};
            snd_outs sndstruct = {};
            fire_struc headstruct = {}, flankstruct = {}, backstruct = {};
            bool activeCrown = false;

            if (args->Simulator == "K")
            {
                // Replica del dinamico: la celda cabeza (para el calculo de pendiente)
                // es el vecino en la direccion de avance; en el borde cae a la misma
                // celda (escenario sin pendiente). perimeterCells = 4*cellSide.
                int r = (cid - 1) / cols, c = (cid - 1) % cols;
                int head_angle = (int)wptr->waz - 90;
                if (head_angle < 0) head_angle += 360;
                head_angle = (int)(std::round(head_angle / 45.0) * 45.0) % 360;
                // vecino en head_angle (convencion de grilla, 0=E,90=N,...)
                int hc = cid;  // fallback: misma celda
                for (int dr = -1; dr <= 1; ++dr)
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        if (dr == 0 && dc == 0) continue;
                        if (neighborAngle(dr, dc) != head_angle) continue;
                        int nr = r + dr, ncl = c + dc;
                        if (nr >= 0 && nr < rows && ncl >= 0 && ncl < cols)
                        {
                            int nid = nr * cols + ncl + 1;
                            if (df[nid - 1].nftype > 0) hc = nid;
                        }
                    }
                calculate_k(&d, &df[hc - 1], (int)perimeterCells, coefs, args, &mainstruct,
                            &sndstruct, &headstruct, &flankstruct, &backstruct, activeCrown, wptr);
            }
            else if (args->Simulator == "C")
                calculate_fbp(&d, coefs, &mainstruct, &sndstruct, &headstruct, &flankstruct,
                              &backstruct, wptr);
            else  // "S" o "P"
                calculate_s(&d, coefs, args, &mainstruct, &sndstruct, &headstruct,
                            &flankstruct, &backstruct, activeCrown, wptr);

            float ros = headstruct.ros * args->HFactor;          // m/min, direccion cabeza
            float inten = mainstruct.sfi;                         // kW/m, Byram superficial
            float fl = mainstruct.fl;                             // m
            // FBP no llena fl en el calculo por celda (depende de la distancia de paso
            // dinamica). Se deriva de la intensidad de Byram con la relacion estandar
            // L = 0.0775 * I^0.46 (Byram), independiente de las tablas de cada kernel.
            if (fl <= 0.0f && inten > 0.0f)
                fl = 0.0775f * std::pow(inten, 0.46f);
            float maxfl = mainstruct.max_flame_length > 0 ? mainstruct.max_flame_length : fl;

            rosL[cid - 1] = ros;
            intenL[cid - 1] = inten;
            flL[cid - 1] = fl;
            crownL[cid - 1] = mainstruct.crown ? 1.0f : 0.0f;
            cfbL[cid - 1] = mainstruct.cfb;
            maxflL[cid - 1] = maxfl;

            // --- ROS DIRECCIONAL por arista (elipse de fuego) ---
            // Para cada vecino combustible j de i, el ROS de i hacia j es el radio de
            // la elipse en la direccion i->j. Usa los mismos semiejes (a,b) y la misma
            // direccion de cabeza que el motor dinamico. El peso de la arista para
            // rutas minimas es el tiempo de transito dist/ROS.
            if (args->StaticGraph && edgeFile.is_open())
            {
                double a = mainstruct.a * args->HFactor;   // semieje mayor (cabeza)
                double b = mainstruct.b * args->FFactor;   // semieje menor (flanco)
                double head = headingCartesian(
                    (args->Simulator == "K" || args->Simulator == "C") ? mainstruct.raz
                                                                       : (double)wptr->waz);
                int r = (cid - 1) / cols, c = (cid - 1) % cols;
                for (int dr = -1; dr <= 1; ++dr)
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        if (dr == 0 && dc == 0) continue;
                        int nr = r + dr, ncl = c + dc;
                        if (nr < 0 || nr >= rows || ncl < 0 || ncl >= cols) continue;
                        int nid = nr * cols + ncl + 1;
                        inputs& dn = df[nid - 1];
                        if (std::strncmp(dn.fueltype, "NF", 2) == 0 ||
                            std::strncmp(dn.fueltype, "ND", 2) == 0 ||
                            dn.fueltype[0] == '\0' || dn.nftype <= 0)
                            continue;  // no se propaga hacia no-combustible
                        int ang = neighborAngle(dr, dc);
                        double offset = ang - head;
                        while (offset < 0) offset += 360.0;
                        while (offset > 360) offset -= 360.0;
                        double rosDir = rhoTheta(offset, a, b) * args->EFactor;  // m/min
                        if (rosDir <= 1e-6) continue;
                        double dist = (dr == 0 || dc == 0) ? cellSide : diagDist;  // m
                        double tmin = dist / rosDir;                               // min
                        // Intensidad de Byram direccional: I propto ROS (Byram clasico),
                        // escalando la intensidad de cabeza por ROS_dir / ROS_cabeza.
                        double fliDir = (ros > 1e-6) ? inten * (rosDir / ros) : 0.0;  // kW/m
                        edgeFile << cid << "," << nid << "," << rosDir << ","
                                 << dist << "," << tmin << "," << fliDir << "\n";
                    }
            }

            // --- Riesgo de breaching por ARCO (misma logica que el motor dinamico) ---
            // Para cada arco bloqueado i->j (barrera entre celdas), el fuego lo cruza si
            // reach = max(k*L, spotFactor*L*U) >= ancho de la barrera. L = longitud de
            // llama de la celda i. Si i puede cruzar hacia algun j, i queda marcada en
            // riesgo de breaching (= salto de barrera posible desde esta celda).
            float breachRisk = 0.0f;
            if (args->BreachFactor > 0.0 && !blockedArcs.empty())
            {
                double L = maxfl > 0.0 ? maxfl : fl;
                double U_ms = wptr->ws / 3.6;
                double reach = std::max((double)args->BreachFactor * L,
                                        (double)args->SpotFactor * L * U_ms);
                int r = (cid - 1) / cols, c = (cid - 1) % cols;
                for (int dr = -1; dr <= 1 && breachRisk == 0.0f; ++dr)
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        if (dr == 0 && dc == 0) continue;
                        int nr = r + dr, ncl = c + dc;
                        if (nr < 0 || nr >= rows || ncl < 0 || ncl >= cols) continue;
                        int nid = nr * cols + ncl + 1;
                        long long key = (long long)cid * nCells + nid;
                        auto it = blockedArcs.find(key);
                        if (it == blockedArcs.end()) continue;  // no hay barrera en este arco
                        if (reach >= it->second)                // la llama supera el ancho
                        {
                            breachRisk = 1.0f;
                            break;
                        }
                    }
            }
            breachL[cid - 1] = breachRisk;
        }

        // Escribir capas de esta hora
        std::string tag = hourTag(weatherData[wp + 1], wp);
        auto W = [&](const std::string& name, std::vector<float>& buf) {
            writeAsc(outDir + name + "_" + tag + ".asc", buf, rows, cols, xllcorner,
                     yllcorner, cellSide);
        };
        W("Intensity", intenL);
        W("ROS", rosL);
        W("FlameLength", flL);
        W("MaxFlameLength", maxflL);
        W("CrownFire", crownL);
        W("CrownFractionBurn", cfbL);
        if (args->BreachFactor > 0.0 && !blockedArcs.empty())
            W("BreachRisk", breachL);

        if (edgeFile.is_open())
        {
            edgeFile.close();
            std::cout << "[sim-static] escenario " << (wp + 1) << "/" << WPeriods << " (" << tag
                      << ") -> capas + grafo direccional (FireGraph_" << tag << ".csv)\n";
        }
        else
        {
            std::cout << "[sim-static] escenario " << (wp + 1) << "/" << WPeriods << " (" << tag
                      << ") -> capas escritas\n";
        }
    }

    std::cout << "[sim-static] listo. Salidas en " << outDir << "\n";
    return 0;
}
