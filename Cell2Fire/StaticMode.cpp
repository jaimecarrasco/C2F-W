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
              int cols, double cellSide, double xllcorner, double yllcorner,
              std::vector<std::vector<std::string>>& weatherData,
              const std::unordered_map<int, double>& barrierWidth,
              const std::unordered_map<std::string, int>& /*fuelLookup*/)
{
    const int WPeriods = static_cast<int>(weatherData.size()) - 1;  // fila 0 = header
    if (WPeriods < 1)
    {
        std::cerr << "[sim-static] Weather.csv sin filas de datos.\n";
        return 1;
    }

    // Parsear TODAS las filas de clima una sola vez (reutiliza el parser del motor).
    CSVReader reader(args->InFolder + "/Weather.csv", ",");
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

    // Loop temporal: una capa-conjunto por fila de clima (hora observada / nowcast /
    // escenario futuro -- cada fila es un mapa estatico independiente).
    for (int wp = 0; wp < WPeriods; ++wp)
    {
        weatherDF* wptr = &wdf[wp];

        for (int cid = 1; cid <= nCells; ++cid)
        {
            inputs& d = df[cid - 1];
            // No-combustible: NODATA en todas las capas
            if (std::strcmp(d.fueltype, NF) == 0 || std::strcmp(d.fueltype, ND) == 0 ||
                d.nftype <= 0)
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
                calculate_k(&d, &d, 8, coefs, args, &mainstruct, &sndstruct, &headstruct,
                            &flankstruct, &backstruct, activeCrown, wptr);
            else if (args->Simulator == "C")
                calculate_fbp(&d, coefs, &mainstruct, &sndstruct, &headstruct, &flankstruct,
                              &backstruct, wptr);
            else  // "S" o "P"
                calculate_s(&d, coefs, args, &mainstruct, &sndstruct, &headstruct,
                            &flankstruct, &backstruct, activeCrown, wptr);

            float ros = headstruct.ros * args->HFactor;          // m/min, direccion cabeza
            float inten = mainstruct.sfi;                         // kW/m, Byram superficial
            float fl = mainstruct.fl;                             // m
            float maxfl = mainstruct.max_flame_length > 0 ? mainstruct.max_flame_length : fl;

            rosL[cid - 1] = ros;
            intenL[cid - 1] = inten;
            flL[cid - 1] = fl;
            crownL[cid - 1] = mainstruct.crown ? 1.0f : 0.0f;
            cfbL[cid - 1] = mainstruct.cfb;
            maxflL[cid - 1] = maxfl;

            // --- Riesgo de breaching por celda ---
            // Celda combustible adyacente a una barrera: cruza si k*FL >= ancho barrera,
            // mismo criterio determinista que el motor dinamico. Alerta -> 1, si no 0.
            float breachRisk = 0.0f;
            if (args->BreachFactor > 0.0 && !barrierWidth.empty())
            {
                int r = (cid - 1) / cols, c = (cid - 1) % cols;
                for (int dr = -1; dr <= 1 && breachRisk == 0.0f; ++dr)
                    for (int dc = -1; dc <= 1; ++dc)
                    {
                        if (dr == 0 && dc == 0)
                            continue;
                        int nr = r + dr, ncl = c + dc;
                        if (nr < 0 || nr >= rows || ncl < 0 || ncl >= cols)
                            continue;
                        int nid = nr * cols + ncl + 1;
                        auto it = barrierWidth.find(nid);
                        if (it == barrierWidth.end())
                            continue;
                        double reach = args->BreachFactor * maxfl;  // alcance de llama (m)
                        if (reach >= it->second)
                        {
                            breachRisk = 1.0f;  // podria saltar esta barrera vecina
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
        if (args->BreachFactor > 0.0 && !barrierWidth.empty())
            W("BreachRisk", breachL);

        std::cout << "[sim-static] escenario " << (wp + 1) << "/" << WPeriods << " (" << tag
                  << ") -> capas escritas\n";
    }

    std::cout << "[sim-static] listo. Salidas en " << outDir << "\n";
    return 0;
}
