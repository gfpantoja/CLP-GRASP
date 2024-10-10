#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <limits>
#include <chrono>
#include <list>
#include <random>
#include <numeric>
#include <omp.h>
#include <cfloat>

using namespace std;

// Variables globales

default_random_engine generator{ 123 };
uniform_real_distribution<double> dist01(0.0, 1.0);
static double const error = 0.0001;
static double errorArea = 100;
static double tolerancia = 0; // tolerancia para que una caja esté dentro de polígono
static auto tIni = chrono::high_resolution_clock::now();
static double duracion = 0;
int holi = 0;
static int ContenedorDimx = 0;
static int ContenedorDimy = 0;
static int ContenedorDimz = 0;
static double ContenedorVol = 0;

static int r_estabilidad = 0; // 0: Sin estabilidad, 1: estabilidad parcial, 2: full support
static int r_juntarEspacios = 0; // 0: Sin juntar ni expandir, 1: juntando sin expandir, 2: juntando y expandiendo
static bool r_maxPresionItems = false; // 0: No hay restricción, 1: si hay restricción
static int r_multidrop = 0; // 0: No hay restricción, 1: Visibilidad, 2: Capas con 0 (Junqueira con d_i = 0)

static int maxNBloques = 500; // Número máximo de bloques que se generan

class Alpha {
public:

    // Parámetros

    double val, funObj, frecuencia;

    // Constructor

    Alpha(double const& p_val) {
        val = p_val;
        funObj = 0;
        frecuencia = 0;
    }

    // Métodos

    void Actualizar(double const& fo) {
        ++frecuencia;
        funObj += (fo - funObj) / frecuencia;
    }

    // Operadores

    bool const operator<(Alpha& otro) {
        return funObj > otro.funObj;
    }
};
vector<Alpha> alphas({ Alpha(0.1), Alpha(0.2), Alpha(0.3), Alpha(0.4), Alpha(0.5), Alpha(0.6), Alpha(0.7), Alpha(0.8), Alpha(0.9), Alpha(1.0) });
vector<double> pesos_acum(alphas.size(), 0);
static double threshold_val = 0.0;
static int threshold_Iter = 0;
const static int threshold_MaxIter = 5;

double totalVolumen = 0;

// Clases

class Point {
public:

    // Parámetros

    int x, y;

    // Constructor

    Point() {
        x = 0; y = 0;
    }
    Point(int const& px, int const& py) {
        x = px;
        y = py;
    }
    Point(Point const& p) {
        x = p.x;
        y = p.y;
    }

    // Métodos

    const double distanciaP(Point const& p) {
        return abs(x - p.x) + abs(y - p.y);
    }

    // Operador

    const bool operator==(Point const& otro) {
        return x == otro.x && y == otro.y;
    }
    const bool operator<(Point const& otro) {
        if (x == otro.x) return y < otro.y;
        return x < otro.x;
    }
};
class Linea {
public:

    // Parámetros

    double A, B, C;

    // Constructor

    Linea() {
        A = 0; B = 0; C = 0;
    }
    Linea(Linea const& l) {
        A = l.A;
        B = l.B;
        C = l.C;
    }
    Linea(Point const& p1, Point const& p2) {
        A = p1.y - p2.y;
        B = p2.x - p1.x;
        double n = sqrtf(A * A + B * B);
        if (n > error) {
            A /= n;
            B /= n;
            C = -p1.x * A - p1.y * B;
        }
        else {
            A = 0.0;
            B = 0.0;
            C = 0.0;
        }
    }
    Linea(double const& x1, double const& y1, double const& x2, double const& y2) {
        A = y1 - y2;
        B = x2 - x1;
        double n = sqrtf(A * A + B * B);
        if (n > error) {
            A /= n;
            B /= n;
            C = -x1 * A - y1 * B;
        }
        else {
            A = 0.0;
            B = 0.0;
            C = 0.0;
        }
    }

    // Métodos

    const double distancia(double const& cx, double const& cy) {
        return A * cx + B * cy + C;
    }
    const double distancia(Point const& p) {
        return A * p.x + B * p.y + C;
    }
    void neg() {
        A = -A; B = -B; C = -C;
    }
};
class Box {
public:

    // Parámetros

    int indBoxType;
    int x, y, z;
    double presion, soporteNecesario;

    // Constructor

    Box() {
        indBoxType = -100;
    }
    Box(Box const& b) {
        indBoxType = b.indBoxType;
        x = b.x;
        y = b.y;
        z = b.z;
        presion = b.presion;
        soporteNecesario = b.soporteNecesario;
    }
    Box(int const& p_indBT, int const& p_x, int const& p_y, int const& p_z) {
        indBoxType = p_indBT;
        x = p_x;
        y = p_y;
        z = p_z;
        presion = 0.0;
        soporteNecesario = 0.0;
    }
    Box(int const& p_indBT, int const& p_x, int const& p_y, int const& p_z, double const& p_presion, double const& masa) {
        indBoxType = p_indBT;
        x = p_x;
        y = p_y;
        z = p_z;
        presion = p_presion;
        soporteNecesario = masa / (double)(p_x * p_y);
    }

    // Operadores

    bool const operator==(Box const& otra) {
        return (x == otra.x && y == otra.y && z == otra.z);
    }
    bool const operator<(Box const& otra) {
        if (x == otra.x) {
            if (y == otra.y) return z < otra.z;
            return y < otra.y;
        }
        return x < otra.x;
    }
};
class BoxType {
public:

    // Parámetros

    int cliente;
    int minDx, minDy, minDz;
    int id, cantidad0, cantidad, volumenReal, idReal;
    double volumen, masa;
    vector<Box> boxes;

private:

    void FullRotation(int const& p_dimx, int const& p_dimy, int const& p_dimz, int const& p_rotx, int const& p_roty, int const& p_rotz) {
        // Rx

        if (p_rotx == 1 || p_rotx == 3) boxes.push_back(Box(id, p_dimz, p_dimy, p_dimx));
        if (p_rotx == 2 || p_rotx == 3) boxes.push_back(Box(id, p_dimy, p_dimz, p_dimx));

        // Ry

        if (p_roty == 1 || p_roty == 3) boxes.push_back(Box(id, p_dimx, p_dimz, p_dimy));
        if (p_roty == 2 || p_roty == 3) boxes.push_back(Box(id, p_dimz, p_dimx, p_dimy));

        // Rz

        if (p_rotz == 1 || p_rotz == 3) boxes.push_back(Box(id, p_dimx, p_dimy, p_dimz));
        if (p_rotz == 2 || p_rotz == 3) boxes.push_back(Box(id, p_dimy, p_dimx, p_dimz));

        // Eliminación de cajas iguales

        sort(boxes.begin(), boxes.end());
        boxes.resize(distance(boxes.begin(), unique(boxes.begin(), boxes.end())));
    }
    void FullRotationPresion(int const& p_dimx, int const& p_dimy, int const& p_dimz, int const& p_rotx, int const& p_roty, int const& p_rotz, double& p_presionx, double& p_presiony, double& p_presionz) {
        // Se ajustan valores de presión

        if (p_dimx == p_dimy)
        {
            if (p_presionx > p_presiony) p_presiony = p_presionx;
            else p_presionx = p_presiony;
        }
        if (p_dimx == p_dimz)
        {
            if (p_presionx > p_presionz) p_presionz = p_presionx;
            else p_presionx = p_presionz;
        }
        if (p_dimy == p_dimz)
        {
            if (p_presiony > p_presionz) p_presionz = p_presiony;
            else p_presiony = p_presionz;
        }

        // Rx

        if (p_rotx == 1 || p_rotx == 3) boxes.push_back(Box(id, p_dimz, p_dimy, p_dimx, p_presionx, masa));
        if (p_rotx == 2 || p_rotx == 3) boxes.push_back(Box(id, p_dimy, p_dimz, p_dimx, p_presionx, masa));


        // Ry

        if (p_roty == 1 || p_roty == 3) boxes.push_back(Box(id, p_dimx, p_dimz, p_dimy, p_presiony, masa));
        if (p_roty == 2 || p_roty == 3) boxes.push_back(Box(id, p_dimz, p_dimx, p_dimy, p_presiony, masa));

        // Rz

        if (p_rotz == 1 || p_rotz == 3) boxes.push_back(Box(id, p_dimx, p_dimy, p_dimz, p_presionz, masa));
        if (p_rotz == 2 || p_rotz == 3) boxes.push_back(Box(id, p_dimy, p_dimx, p_dimz, p_presionz, masa));

        // Eliminación de cajas iguales

        sort(boxes.begin(), boxes.end());
        boxes.resize(distance(boxes.begin(), unique(boxes.begin(), boxes.end())));
    }
    void DetermineDmins() {
        minDx = boxes.front().x;
        minDy = boxes.front().y;
        minDz = boxes.front().z;
        for (vector<Box>::iterator b_it = boxes.begin() + 1; b_it < boxes.end(); ++b_it) {
            if ((*b_it).x < minDx) minDx = (*b_it).x;
            if ((*b_it).y < minDy) minDy = (*b_it).y;
            if ((*b_it).z < minDz) minDz = (*b_it).z;
        }
    }

public:

    // Constructores

    BoxType() {
        boxes = vector<Box>(0);
    }
    BoxType(BoxType const& bt) {
        cliente = bt.cliente;
        minDx = bt.minDx;
        minDy = bt.minDy;
        minDz = bt.minDz;
        id = bt.id;
        idReal = bt.idReal;
        cantidad = bt.cantidad;
        cantidad0 = bt.cantidad0;
        volumen = bt.volumen;
        volumenReal = bt.volumenReal;
        masa = bt.masa;
        boxes = bt.boxes;
    }
    BoxType(int const& p_id, int const& p_cantidad, int const& p_dimx, int const& p_dimy, int const& p_dimz, int const& p_rotx, int const& p_roty, int const& p_rotz, double const& p_masa = 0.0f, int const& p_cliente = 0) {
        id = p_id;
        idReal = p_id;
        cantidad = p_cantidad;
        cantidad0 = p_cantidad;
        cliente = p_cliente;
        masa = p_masa;
        volumenReal = p_dimx * p_dimy * p_dimz;
        volumen = (double)volumenReal / ContenedorVol * 100.0;
        boxes.reserve(6);
        FullRotation(p_dimx, p_dimy, p_dimz, p_rotx, p_roty, p_rotz);
        DetermineDmins();
    }
    BoxType(int const& p_id, int const& p_cantidad, int const& p_dimx, int const& p_dimy, int const& p_dimz, int const& p_rotx, int const& p_roty, int const& p_rotz, double const& p_masa, double& p_presionx, double& p_presiony, double& p_presionz, int const& p_cliente = 0) {
        id = p_id;
        idReal = p_id;
        cantidad = p_cantidad;
        cantidad0 = p_cantidad;
        cliente = p_cliente;
        masa = p_masa;
        volumenReal = p_dimx * p_dimy * p_dimz;
        volumen = (double)volumenReal / ContenedorVol * 100.0;
        boxes.reserve(6);
        FullRotationPresion(p_dimx, p_dimy, p_dimz, p_rotx, p_roty, p_rotz, p_presionx, p_presiony, p_presionz);
        DetermineDmins();
    }

    // Métodos

    void ActualizarIndice(int const& p_id) {
        id = p_id;
        for (vector<Box>::iterator b_it = boxes.begin(); b_it < boxes.end(); ++b_it) (*b_it).indBoxType = p_id;
    }

    // Operadores

    bool const operator<(BoxType const& otro) {
        return id < otro.id;
    }
    bool const operator==(BoxType const& otro) {
        return id == otro.id;
    }
};
class MaximalSpace {
public:

    // Parámetros

    int indEsquina; // 1: (-x, -y, -z), 2: (+x, -y, -z), 3: (-x +y, -z), 4: (+x, +y, -z) -> los demás con +z
    int metrica, sumDim;
    int x1, y1, z1, x2, y2, z2, dx, dy, dz;

    bool cambio; // Indica si hay cambio o no para determinar los indices de las cajas
    vector<int> indCajas; // Indices de las cajas empacadas sobre la que se forma el EM
    double maxSoporte; // máximo soporte de las cajas debajo del EM

    // Constructor

    MaximalSpace() {
        x1 = 0;
        y1 = 0;
        z1 = 0;
        y2 = ContenedorDimy;
        x2 = ContenedorDimx;
        z2 = ContenedorDimz;
        dx = ContenedorDimx;
        dy = ContenedorDimy;
        dz = ContenedorDimz;
        indEsquina = 1;
        indCajas = vector<int>(0);
        maxSoporte = FLT_MAX;
    }
    MaximalSpace(int const& p_e1x, int const& p_e1y, int const& p_e1z, int const& p_e2x, int const& p_e2y, int const& p_e2z) {
        x1 = p_e1x;
        y1 = p_e1y;
        z1 = p_e1z;
        x2 = p_e2x;
        y2 = p_e2y;
        z2 = p_e2z;
        dx = p_e2x - p_e1x;
        dy = p_e2y - p_e1y;
        dz = p_e2z - p_e1z;
        sumDim = dx + dy + dz;
        indCajas = vector<int>(0);
        maxSoporte = FLT_MAX;
        cambio = true;
        DeterminarEsquina();
    }

    // Metodo

    void DeterminarEsquina() {
        /*
        metrica = x1 + y1 + z1;
        indEsquina = 1;
        if (metrica > 0) {
            int miy = ContenedorDimy - y2;
            int val = x1 + miy + z1;
            if (metrica > val) {
                metrica = val;
                indEsquina = 3;
                if (metrica == 0) return;
            }
        }
        */
        if (r_multidrop > 0) {
            metrica = x1 + y1 + z1;
            indEsquina = 1;
            /*
            if (r_estabilidad == 0) { // Se consideran las 4 esquinas hacia la izquierda (-x)
                metrica = x1 + y1 + z1;
                indEsquina = 1;
                if (metrica > 0) {
                    int miy = ContenedorDimy - y2;
                    int val = x1 + miy + z1;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 3;
                        if (metrica == 0) return;
                    }
                    int miz = ContenedorDimz - z2;
                    val = x1 + y1 + miz;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 5;
                        if (metrica == 0) return;
                    }
                    val = x1 + miy + miz;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 7;
                        if (metrica == 0) return;
                    }
                }
            }
            else { // Se consideran las 2 esquinas hacia la izquierda y hacia abajo (-x, -z)
                metrica = x1 + y1 + z1;
                indEsquina = 1;
                if (metrica > 0) {
                    int miy = ContenedorDimy - y2;
                    int val = x1 + miy + z1;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 3;
                        if (metrica == 0) return;
                    }
                }
            }
            */
        }
        else {
            if (r_estabilidad == 0) { // Se consideran las 8 esquinas
                metrica = x1 + y1 + z1;
                indEsquina = 1;
                if (metrica > 0) {
                    int mix = ContenedorDimx - x2;
                    int val = mix + y1 + z1;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 2;
                        if (metrica == 0) return;
                    }
                    int miy = ContenedorDimy - y2;
                    val = x1 + miy + z1;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 3;
                        if (metrica == 0) return;
                    }
                    val = mix + miy + z1;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 4;
                        if (metrica == 0) return;
                    }
                    int miz = ContenedorDimz - z2;
                    val = x1 + y1 + miz;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 5;
                        if (metrica == 0) return;
                    }
                    val = mix + y1 + miz;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 6;
                        if (metrica == 0) return;
                    }
                    val = x1 + miy + miz;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 7;
                        if (metrica == 0) return;
                    }
                    val = mix + miy + miz;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 8;
                        if (metrica == 0) return;
                    }
                }
            }
            else { // Se consideran las 4 esquinas hacia abajo (-z)
                metrica = x1 + y1 + z1;
                indEsquina = 1;
                if (metrica > 0) {
                    int mix = ContenedorDimx - x2;
                    int val = mix + y1 + z1;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 2;
                        if (metrica == 0) return;
                    }
                    int miy = ContenedorDimy - y2;
                    val = x1 + miy + z1;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 3;
                        if (metrica == 0) return;
                    }
                    val = mix + miy + z1;
                    if (metrica > val) {
                        metrica = val;
                        indEsquina = 4;
                        if (metrica == 0) return;
                    }
                }
            }
        }

    }

    // Operadores

    bool const operator<(MaximalSpace const& otro) {
        if (x1 == otro.x1) {
            if (z1 == otro.z1) {
                return metrica < otro.metrica;
            }
            return z1 > otro.z1;
        }
        return x1 < otro.x1;
    }
};
class BloqueAux {
public:

    // Parámetros

    int x1, y1, z1;

    // Constructor

    BloqueAux(int const& p_x1, int const& p_y1, int const& p_z1) {
        x1 = p_x1;
        y1 = p_y1;
        z1 = p_z1;
    }
};
class Bloque {
public:

    // Parámetros

    int x1, y1, z1, x2, y2, z2, dx, dy, dz, esquina, bf1, bf2, bf3, bfx, bfy, bfz, vol;
    BoxType bt;
    Box b;
    vector<BloqueAux> cajas;
    int btq;
    vector<int> indPB_EM; // indices de las posiciones de las cajas empacadas debajo del EM

    // Constructor

    Bloque() {

    }
    Bloque(BoxType const& p_bt, Box const& p_b, MaximalSpace const& p_e) {
        dx = p_b.x;
        dy = p_b.y;
        dz = p_b.z;
        vol = p_bt.volumenReal;
        indPB_EM = p_e.indCajas;
        esquina = p_e.indEsquina;
        if (p_e.indEsquina == 1) {
            x1 = p_e.x1;
            y1 = p_e.y1;
            z1 = p_e.z1;
            x2 = p_e.x1 + p_b.x;
            y2 = p_e.y1 + p_b.y;
            z2 = p_e.z1 + p_b.z;
        }
        else if (p_e.indEsquina == 2) {
            x1 = p_e.x2 - p_b.x;
            y1 = p_e.y1;
            z1 = p_e.z1;
            x2 = p_e.x2;
            y2 = p_e.y1 + p_b.y;
            z2 = p_e.z1 + p_b.z;
        }
        else if (p_e.indEsquina == 3) {
            x1 = p_e.x1;
            y1 = p_e.y2 - p_b.y;
            z1 = p_e.z1;
            x2 = p_e.x1 + p_b.x;
            y2 = p_e.y2;
            z2 = p_e.z1 + p_b.z;
        }
        else if (p_e.indEsquina == 4) {
            x1 = p_e.x2 - p_b.x;
            y1 = p_e.y2 - p_b.y;
            z1 = p_e.z1;
            x2 = p_e.x2;
            y2 = p_e.y2;
            z2 = p_e.z1 + p_b.z;
        }
        else if (p_e.indEsquina == 5) {
            x1 = p_e.x1;
            y1 = p_e.y1;
            z1 = p_e.z2 - p_b.z;
            x2 = p_e.x1 + p_b.x;
            y2 = p_e.y1 + p_b.y;
            z2 = p_e.z2;
        }
        else if (p_e.indEsquina == 6) {
            x1 = p_e.x2 - p_b.x;
            y1 = p_e.y1;
            z1 = p_e.z2 - p_b.z;
            x2 = p_e.x2;
            y2 = p_e.y1 + p_b.y;
            z2 = p_e.z2;
        }
        else if (p_e.indEsquina == 7) {
            x1 = p_e.x1;
            y1 = p_e.y2 - p_b.y;
            z1 = p_e.z2 - p_b.z;
            x2 = p_e.x1 + p_b.x;
            y2 = p_e.y2;
            z2 = p_e.z2;
        }
        else { // indEsquina = 8
            x1 = p_e.x2 - p_b.x;
            y1 = p_e.y2 - p_b.y;
            z1 = p_e.z2 - p_b.z;
            x2 = p_e.x2;
            y2 = p_e.y2;
            z2 = p_e.z2;
        }
        b = p_b;
        cajas = vector<BloqueAux>(1, BloqueAux(x1, y1, z1));
        bt = p_bt;
        btq = 1;
        bfx = p_e.dx - dx;
        bf1 = bfx;
        bfy = p_e.dy - dy;
        bf2 = bfy;
        bfz = p_e.dz - dz;
        bf3 = bfz;
        if (bf1 > bf2) {
            int temp = bf1;
            bf1 = bf2;
            bf2 = temp;
        }
        if (bf1 > bf3) {
            int temp = bf1;
            bf1 = bf3;
            bf3 = temp;
        }
        if (bf2 > bf3) {
            int temp = bf2;
            bf2 = bf3;
            bf3 = temp;
        }
    }
    Bloque(BoxType const& p_bt, Box const& p_b, MaximalSpace const& p_e, int const& nx, int const& ny, int const& nz) {
        dx = p_b.x * nx;
        dy = p_b.y * ny;
        dz = p_b.z * nz;
        vol = p_bt.volumenReal * nx * ny * nz;
        indPB_EM = p_e.indCajas;
        esquina = p_e.indEsquina;
        if (p_e.indEsquina == 1) {
            x1 = p_e.x1;
            y1 = p_e.y1;
            z1 = p_e.z1;
            x2 = p_e.x1 + p_b.x * nx;
            y2 = p_e.y1 + p_b.y * ny;
            z2 = p_e.z1 + p_b.z * nz;
        }
        else if (p_e.indEsquina == 2) {
            x1 = p_e.x2 - p_b.x * nx;
            y1 = p_e.y1;
            z1 = p_e.z1;
            x2 = p_e.x2;
            y2 = p_e.y1 + p_b.y * ny;
            z2 = p_e.z1 + p_b.z * nz;
        }
        else if (p_e.indEsquina == 3) {
            x1 = p_e.x1;
            y1 = p_e.y2 - p_b.y * ny;
            z1 = p_e.z1;
            x2 = p_e.x1 + p_b.x * nx;
            y2 = p_e.y2;
            z2 = p_e.z1 + p_b.z * nz;
        }
        else if (p_e.indEsquina == 4) {
            x1 = p_e.x2 - p_b.x * nx;
            y1 = p_e.y2 - p_b.y * ny;
            z1 = p_e.z1;
            x2 = p_e.x2;
            y2 = p_e.y2;
            z2 = p_e.z1 + p_b.z * nz;
        }
        else if (p_e.indEsquina == 5) {
            x1 = p_e.x1;
            y1 = p_e.y1;
            z1 = p_e.z2 - p_b.z * nz;
            x2 = p_e.x1 + p_b.x * nx;
            y2 = p_e.y1 + p_b.y * ny;
            z2 = p_e.z2;
        }
        else if (p_e.indEsquina == 6) {
            x1 = p_e.x2 - p_b.x * nx;
            y1 = p_e.y1;
            z1 = p_e.z2 - p_b.z * nz;
            x2 = p_e.x2;
            y2 = p_e.y1 + p_b.y * ny;
            z2 = p_e.z2;
        }
        else if (p_e.indEsquina == 7) {
            x1 = p_e.x1;
            y1 = p_e.y2 - p_b.y * ny;
            z1 = p_e.z2 - p_b.z * nz;
            x2 = p_e.x1 + p_b.x * nx;
            y2 = p_e.y2;
            z2 = p_e.z2;
        }
        else { // indEsquina = 8
            x1 = p_e.x2 - p_b.x * nx;
            y1 = p_e.y2 - p_b.y * ny;
            z1 = p_e.z2 - p_b.z * nz;
            x2 = p_e.x2;
            y2 = p_e.y2;
            z2 = p_e.z2;
        }
        bfx = p_e.dx - dx;
        bf1 = bfx;
        bfy = p_e.dy - dy;
        bf2 = bfy;
        bfz = p_e.dz - dz;
        bf3 = bfz;
        if (bf1 > bf2) {
            int temp = bf1;
            bf1 = bf2;
            bf2 = temp;
        }
        if (bf1 > bf3) {
            int temp = bf1;
            bf1 = bf3;
            bf3 = temp;
        }
        if (bf2 > bf3) {
            int temp = bf2;
            bf2 = bf3;
            bf3 = temp;
        }
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    cajas.push_back(BloqueAux(x1 + p_b.x * x, y1 + p_b.y * y, z1 + p_b.z * z));
                }
            }
        }
        b = p_b;
        bt = p_bt;
        btq = nx * ny * nz;
    }

    // Operator

    bool const operator==(Bloque const& otro) {
        return dx == otro.dx && dy == otro.dz && b == otro.b;
    }
    bool const operator<(Bloque const& otro) {
        if (dx == otro.dx) {
            if (dy == otro.dy) {
                if (dz == otro.dz) {
                    return b < otro.b;
                }
                return dz < otro.dz;
            }
            return dy < otro.dy;
        }
        return dx < otro.dx;
    }
};
class Soporte {
public:

    // Parámetros

    int indPB; // índice de la caja empacada en la lista del contenedor
    double area; // area de contacto
    double centrox, centroy; // centro del area de contacto
    double masaSoporte; // Esta cambia al agregar cajas encima
    double deltaMasa; // Cambio de masa en el soporte
    double masaMax; // masa que máximo soporta según el area y la presión 

    // Constructor

    Soporte() {

    }
    Soporte(Soporte const& s) {
        indPB = s.indPB;
        area = s.area;
        centrox = s.centrox;
        centroy = s.centroy;
        masaSoporte = s.masaSoporte;
        deltaMasa = s.deltaMasa;
        masaMax = s.masaMax;
    }
    Soporte(int const& p_indPB, int const& p_area, double const& p_cx, double const& p_cy, double const& p_presion) {
        indPB = p_indPB;
        area = (double)p_area;
        centrox = p_cx;
        centroy = p_cy;
        masaMax = p_presion * area;
        masaSoporte = 0;
        deltaMasa = 0;
    }

    // Métodos

    const bool SiSeSoporta() {
        return masaMax >= masaSoporte;
    }
    void ActualizarSoporte(double const& new_Masa) {
        deltaMasa = new_Masa - masaSoporte;
        masaSoporte = new_Masa;
    }
};
class PolAux {

    // Parámetros

private:

    Linea linea;

public:

    vector<Point> grupo;
    Point v1, v2;

    // Constructor

    PolAux(vector<Point> const& p_grupo, Linea const& p_li, Point const& p_v1, Point const& p_v2) {
        grupo = p_grupo;
        linea = p_li;
        v1 = p_v1;
        v2 = p_v2;
    }

    // Métodos

    const int DeterminarAlejado() {
        double dAlejado = linea.distancia(grupo.front());
        int id = 0;
        int i = 1;
        for (vector<Point>::iterator g_it = grupo.begin() + 1; g_it < grupo.end(); ++g_it, ++i) {
            double d = linea.distancia(*g_it);
            if (d < dAlejado) {
                dAlejado = d;
                id = i;
            }
        }
        return id;
    }
};
class Poligono {
public:

    // Parámetros

    double xmin, xmax, ymin, ymax;
    vector<Linea> lineas;

    // Constructor

    Poligono() {
        xmin = 0;
        xmax = 0;
        ymin = 0;
        ymax = 0;
        lineas = vector<Linea>(0);
    }
    Poligono(Poligono const& p) {
        xmin = p.xmin;
        xmax = p.xmax;
        ymin = p.ymin;
        ymax = p.ymax;
        lineas = p.lineas;
    }
    Poligono(vector<Point>& points)
    {
        // Obtener únicos

        sort(points.begin(), points.end());
        points.resize(distance(points.begin(), unique(points.begin(), points.end())));

        // Encontrar 2 vértices extremos

        vector<Point> ch({ points.front(), points.back() });
        xmin = points.front().x;
        xmax = points.back().x;
        Linea l12(points.front(), points.back());
        //points.pop_back();
        //points.erase(points.begin());

        // Se dividen los puntos

        vector<Point> grupo1; grupo1.reserve(points.size());
        vector<Point> grupo2; grupo2.reserve(points.size());
        for (vector<Point>::iterator p_it = points.begin() + 1; p_it < points.end() - 1; ++p_it) {
            double d = l12.distancia(*p_it);
            if (d < -error) grupo2.push_back(*p_it);
            else if (d > error) grupo1.push_back(*p_it);
        }

        // Se crea la lista de PolAux

        vector<PolAux> lista; lista.reserve(points.size());
        if (grupo2.size() > 0) lista.push_back(PolAux(grupo2, l12, ch.front(), ch.back()));
        if (grupo1.size() > 0) {
            l12.neg();
            lista.push_back(PolAux(grupo1, l12, ch.back(), ch.front()));
        }

        while (lista.size() > 0) {
            vector<PolAux> copia = lista;
            lista.clear();
            for (vector<PolAux>::iterator p = copia.begin(); p < copia.end(); ++p) {
                int ind = (*p).DeterminarAlejado();
                ch.push_back((*p).grupo[ind]);
                //(*p).grupo.erase((*p).grupo.begin() + ind);

                // Linea 1

                vector<Point> grupo; grupo.reserve((*p).grupo.size());
                Linea linea((*p).v1, ch.back());
                for (vector<Point>::iterator po = (*p).grupo.begin(); po < (*p).grupo.end(); ++po) {
                    if (linea.distancia(*po) < -error) grupo.push_back(*po);
                }
                if (grupo.size() > 0) lista.push_back(PolAux(grupo, linea, (*p).v1, ch.back()));

                // Linea 2

                grupo.clear();
                linea = Linea(ch.back(), (*p).v2);
                for (vector<Point>::iterator po = (*p).grupo.begin(); po < (*p).grupo.end(); ++po) {
                    if (linea.distancia(*po) < -error) grupo.push_back(*po);
                }
                if (grupo.size() > 0) lista.push_back(PolAux(grupo, linea, ch.back(), (*p).v2));
            }
        }

        // min y max

        ymin = ch.front().y;
        ymax = ch.front().y;
        for (vector<Point>::iterator p = ch.begin() + 1; p < ch.end(); ++p) {
            if (ymin > (*p).y) ymin = (*p).y;
            else if (ymax < (*p).y) ymax = (*p).y;
        }

        // Lineas

        lineas.reserve(ch.size());
        for (int i = 0; i < (int)ch.size() - 1; ++i) {
            Point point1 = ch[i];
            for (int j = i + 1; j < ch.size(); ++j) {
                Point point2 = ch[j];
                double maxD = point1.distanciaP(point2);
                Linea l(point1, point2);
                bool sirve = true;
                for (int k = (int)ch.size() - 1; k >= 0; --k) {
                    if (k != i && k != j) {
                        Point point3 = ch[k];
                        double d = l.distancia(point3);
                        if (d < -error) {
                            sirve = false;
                            break;
                        }
                        if (k > i) {
                            if (abs(d) < error) { // colineal
                                d = point1.distanciaP(point3);
                                if (maxD < d) {
                                    maxD = d;
                                    ch[j] = ch[k];
                                    point2 = ch[j];
                                    l = Linea(point1, point2);
                                }
                                ch.erase(ch.begin() + k);
                                if (k < j) --j;
                            }
                        }
                    }
                }
                if (sirve) {
                    lineas.push_back(l);
                    Point temp = ch[j];
                    ch[j] = ch[i + 1];
                    ch[i + 1] = temp;
                    break;
                }
            }
        }
        lineas.push_back(Linea(ch.back(), ch.front()));
    }

    // Métodos

    const bool puntoDentroPoligono(double const& cx, double const& cy) {
        if (xmin <= cx && cx <= xmax) {
            if (ymin <= cy && cy <= ymax) {
                for (vector<Linea>::iterator l = lineas.begin(); l < lineas.end(); ++l) {
                    if ((*l).distancia(cx, cy) < tolerancia) return false;
                }
                return true;
            }
            return false;
        }
        return false;
    }
};
class PackedBox {
public:

    // Parámetros

    int x1, y1, z1, x2, y2, z2;
    int id, ind, cliente; // el id es el id del boxtype
    double volumen, masa, presion;

    double masaTotal; // Sumatoria de las masas
    double centroMasax, centroMasay;
    vector<Soporte> soportes; // Cajas debajo de esta
    Poligono poligono;

    // Constructor

    PackedBox() {
        soportes = vector<Soporte>(0);
    }
    PackedBox(PackedBox const& p) {
        x1 = p.x1;
        y1 = p.y1;
        z1 = p.z1;
        x2 = p.x2;
        y2 = p.y2;
        z2 = p.z2;
        id = p.id;
        ind = p.ind;
        cliente = p.cliente;
        volumen = p.volumen;
        masa = p.masa;
        presion = p.presion;
        masaTotal = p.masaTotal;
        centroMasax = p.centroMasax;
        centroMasay = p.centroMasay;
        poligono = p.poligono;
        soportes = p.soportes;
    }
    PackedBox(int const& p_ind, Bloque& bloque, BloqueAux const& bloqueAux) { // Se usa para empacar cuando no hay presion
        Inicializar(p_ind, bloque, bloqueAux);
    }
    PackedBox(int const& p_ind, Bloque& bloque, BloqueAux const& bloqueAux, double const& cmx, double const& cmy, Poligono const& p_pol, vector<Soporte> const& p_soportes) {
        Inicializar(p_ind, bloque, bloqueAux);
        masaTotal = masa;
        centroMasax = cmx;
        centroMasay = cmy;
        poligono = p_pol;
        soportes = p_soportes;
    }

    // Metodos

    void ActualizarCentroDeMasa(Soporte const& s) {
        double newMasaTotal = s.deltaMasa + masaTotal;
        centroMasax = (masaTotal * centroMasax + s.deltaMasa * s.centrox) / newMasaTotal;
        centroMasay = (masaTotal * centroMasay + s.deltaMasa * s.centroy) / newMasaTotal;
        masaTotal = newMasaTotal;
    }
    bool ActualizarSoportes() {
        vector<double> resp;
        ObtenerEcuacionesSoporte(resp);
        int i = 0;
        for (vector<Soporte>::iterator s = soportes.begin(); s < soportes.end(); ++s, ++i) {
            (*s).ActualizarSoporte(resp[i]);

            if (!(*s).SiSeSoporta()) {
                return false;
            }
        }
        return true;
    }
    bool centroMasaEnPoligono() {
        return poligono.puntoDentroPoligono(centroMasax, centroMasay);
    }
private:

    void Inicializar(int const& p_ind, Bloque& bloque, BloqueAux const& bloqueAux) {
        ind = p_ind;
        id = bloque.bt.id;
        cliente = bloque.bt.cliente;
        volumen = bloque.bt.volumen;
        masa = bloque.bt.masa;
        presion = bloque.b.presion;
        x1 = bloqueAux.x1;
        y1 = bloqueAux.y1;
        z1 = bloqueAux.z1;
        x2 = x1 + bloque.b.x;
        y2 = y1 + bloque.b.y;
        z2 = z1 + bloque.b.z;
    }
    void ObtenerEcuacionesSoporte(vector<double>& resp) {
        vector<vector<double>> A;
        vector<double> A1; A1.reserve(soportes.size() + 1);
        vector<double> A2; A2.reserve(soportes.size() + 1);
        for (vector<Soporte>::iterator s = soportes.begin(); s < soportes.end(); ++s) {
            A1.push_back(centroMasax - (*s).centrox);
            A2.push_back(centroMasay - (*s).centroy);
        }
        A1.push_back(0);
        A2.push_back(0);
        A.push_back(A1);
        A.push_back(A2);
        vector<double> A3; A3.reserve(soportes.size() + 1);
        for (vector<Soporte>::iterator s1 = soportes.begin(); s1 < soportes.end(); ++s1) {
            A1.clear();
            A2.clear();
            for (vector<Soporte>::iterator s2 = soportes.begin(); s2 < soportes.end(); ++s2) {
                A1.push_back((*s1).centrox - (*s2).centrox);
                A2.push_back((*s1).centroy - (*s2).centroy);
            }
            A1.push_back(((*s1).centrox - centroMasax) * masaTotal);
            A2.push_back(((*s1).centroy - centroMasay) * masaTotal);
            A.push_back(A1);
            A.push_back(A2);
            A3.push_back(1);
        }
        A3.push_back(masaTotal);
        A.push_back(A3); // Sumatoria de fuerzas
        resp.clear(); resp.reserve(soportes.size());
        GaussianElimination(A, resp);
        if (A.size() < (int)A.front().size() - 1) { // Se supone presion igual
            int nColumnas = A.front().size();

            // Se añaden todas las ecuaciones de soporte según área

            int i = 0;
            for (vector<Soporte>::iterator s_it1 = soportes.begin(); s_it1 < soportes.end() - 1; ++s_it1, ++i) {
                int j = i + 1;
                for (vector<Soporte>::iterator s_it2 = s_it1 + 1; s_it2 < soportes.end(); ++s_it2, ++j) {
                    A.push_back(vector<double>(nColumnas, 0));
                    A.back()[i] = 1.0f / (*s_it1).area;
                    A.back()[j] = -1.0f / (*s_it2).area;
                }
            }

            // Eliminación gaussiana 2

            GaussianElimination2(A, resp);
        }
    }
    void GaussianElimination(vector<vector<double>>& A, vector<double>& resp) {

        // La matriz A son los indices de cada variable y el resultado negativo

        int nF = (int)A.front().size() - 1;
        resp.resize(nF);
        int nR = A.size();

        // Obtener matriz triangular

        for (int j = 0; j < nF; ++j) {

            // Encontrar el máximo valor absoluto de la columna hacia abajo

            int imax = j;
            double valmax = abs(A[imax][j]);
            for (int i = j + 1; i < nR; ++i) {
                if (abs(A[i][j]) > valmax) {
                    imax = i;
                    valmax = abs(A[i][j]);
                }
            }
            if (valmax < errorArea) {
                A.erase(A.begin() + j, A.begin() + (int)A.size() - j);

                // Poner 0s

                int i = 1;
                for (vector<vector<double>>::iterator A2 = A.begin() + 1; A2 < A.end(); ++A2, ++i) {
                    vector<double>::iterator A1 = (*A2).begin();
                    for (int j2 = 0; j2 < i; ++j2, ++A1) {
                        (*A1) = 0;
                    }
                }
                return;
            }
            valmax = A[imax][j];

            // Intercambiar fila actual con la de máximo valor

            if (imax != j) {
                vector<double> temp(A[imax]);
                A[imax] = A[j];
                A[j] = temp;
            }

            // Dividir las siguientes filas y columnas según un factor para tener una matriz triangular

            int i2 = j + 1;
            for (int i = j + 1; i < nR; ++i, ++i2) {
                double factor = A[i][j] / valmax;
                if (abs(factor) < errorArea) continue;
                int j2 = j + 1;
                for (int k = j + 1; k < nF + 1; ++k, ++j2) {
                    A[i][k] -= A[j][j2] * factor;
                }
                A[i2][j] = 0.0;
            }

            // Poner el pivote en 1

            A[j][j] = 1;
            for (int k = j + 1; k < nF + 1; ++k) {
                A[j][k] /= valmax;
            }
        }
        A.resize(nF);

        // Obtener respuesta a partir de resolver las ecuaciones con la matriz triangular de abajo para arriba

        int ii = nF - 1;
        for (int k = (int)A.size() - 1; k >= 0; --k, --ii) {
            resp[ii] = A[k][nF];
            int j = ii + 1;
            for (int m = ii + 1; m < nF; ++m, ++j) {
                resp[ii] -= A[k][m] * resp[j];
            }
        }
    }
    void GaussianElimination2(vector<vector<double>>& A, vector<double>& resp) {

        // La matriz A son los indices de cada variable y el resultado negativo

        int nF = (int)A.front().size() - 1;// Número de variables
        int nR = A.size();// Número de ecuaciones 

        // Obtener matriz triangular

        for (int j = 0; j < nF; ++j)
        {
            // Si el número actual es diferente de cero se deja ese, sino se determina el máximo valor absoluto. Esto para ir en orden de las ecuaciones

            int imax = j;
            double valmax = abs(A[imax][j]);
            if (valmax < errorArea) {
                for (int i = j + 1; i < nR; ++i) {
                    if (abs(A[i][j]) > valmax) {
                        imax = i;
                        valmax = abs(A[i][j]);
                    }
                }
            }
            if (valmax < errorArea) {
                // MACHETAZO
                for (vector<double>::iterator r_it = resp.begin(); r_it < resp.end(); ++r_it) {
                    (*r_it) = numeric_limits<double>::max();
                }
                cout << "Matriz singular" << endl;
                return;
            }
            valmax = A[imax][j];

            // Intercambiar fila actual con la de máximo valor

            if (imax != j)
            {
                vector<double> temp(A[imax]);
                A[imax] = A[j];
                A[j] = temp;
            }

            // Dividir las siguientes filas y columnas según un factor para tener una matriz triangular

            int i2 = j + 1;
            for (int i = j + 1; i < nR; ++i)
            {
                double factor = A[i][j] / valmax;
                if (abs(factor) < errorArea) continue;
                int j2 = j + 1;
                for (int k = j + 1; k < nF + 1; ++k)
                {
                    A[i][k] -= A[j][j2] * factor;
                    ++j2;
                }
                A[i2][j] = 0.0;
                ++i2;
            }

            // Poner el pivote en 1

            A[j][j] = 1;
            for (int k = j + 1; k < nF + 1; ++k)
            {
                A[j][k] /= valmax;
            }
        }
        A.resize(nF);

        // Obtener respuesta a partir de resolver las ecuaciones con la matriz triangular de abajo para arriba

        int ii = nF - 1;
        for (int k = (int)A.size() - 1; k >= 0; --k, --ii) {
            resp[ii] = A[k][nF];
            int j = ii + 1;
            for (int m = ii + 1; m < nF; ++m, ++j) {
                resp[ii] -= A[k][m] * resp[j];
            }
        }
    }
};
bool const Orden_Fondo_Altura(MaximalSpace const& m1, MaximalSpace const& m2) {
    if (m1.x1 == m2.x1) {
        if (m1.z1 == m2.z1) {
            return m1.metrica < m2.metrica;
        }
        return m1.z1 > m2.z1;
    }
    return m1.x1 < m2.x1;
}
bool const Orden_Fondo_Bajo(MaximalSpace const& m1, MaximalSpace const& m2) {
    if (m1.x1 == m2.x1) {
        if (m1.z1 == m2.z1) {
            return m1.metrica < m2.metrica;
        }
        return m1.z1 < m2.z1;
    }
    return m1.x1 < m2.x1;
}
bool const Orden_Altura_Fondo(MaximalSpace const& m1, MaximalSpace const& m2) {
    if (m1.z1 == m2.z1) {
        if (m1.x1 == m2.x1) {
            return m1.metrica < m2.metrica;
        }
        return m1.x1 < m2.x1;
    }
    return m1.z1 > m2.z1;
}
bool const Orden_BestFit_MaxVol(Bloque const& b1, Bloque const& b2) {
    /*
    if (b1.bf1 == b2.bf1) {
        if (b1.bf2 == b2.bf2) {
            if (b1.bf3 == b2.bf3) {
                return b1.vol > b2.vol;
            }
            return b1.bf3 < b2.bf3;
        }
        return b1.bf2 < b2.bf2;
    }
    */
    return b1.bf1 < b2.bf1;
}
bool const Orden_MaxVol_BestFit(Bloque const& b1, Bloque const& b2) {
    /*
    if (b1.vol == b2.vol) {
        if (b1.bf1 == b2.bf1) {
            if (b1.bf2 == b2.bf2) {
                return b1.bf3 < b2.bf3;
            }
            return b1.bf2 < b2.bf2;
        }
        return b1.bf1 < b2.bf1;
    }
    */
    return b1.vol > b2.vol;
}
class Container {

    // Parámetros

public:
    int minDx, minDy, minDz;
    int clienteActual;
    double alpha, alphaEM;
    vector<MaximalSpace> espacios;
    double utilizacion, volEmpacado;
    vector<BoxType> boxest;
    vector<PackedBox> empacados;
    bool hayCajasPorEmpacar;

    // Constructor

    Container() {
        espacios = vector<MaximalSpace>(0);
        boxest = vector<BoxType>(0);
        empacados = vector<PackedBox>(0);
    }
    Container(vector<BoxType> const& p_boxest, double const& p_alpha, double const& p_alphaEM = 0) {
        boxest = p_boxest;
        alpha = p_alpha;
        alphaEM = p_alphaEM;
        ActualizarMinD();
        utilizacion = 0.0;
        espacios = vector<MaximalSpace>(1);
        clienteActual = 0;
        hayCajasPorEmpacar = true;
    }
    Container(Container& c, bool const& copiaCompleta) {
        if (copiaCompleta) {
            alpha = c.alpha;
            alpha = c.alphaEM;
            minDx = c.minDx;
            minDy = c.minDy;
            minDz = c.minDz;
            boxest = c.boxest;
            utilizacion = 0.0;
            espacios = vector<MaximalSpace>(1);
            empacados.clear();
            clienteActual = 0;
            hayCajasPorEmpacar = true;
        }
        else {
            utilizacion = c.utilizacion;
            empacados = c.empacados;
            //boxest = c.boxest;
            hayCajasPorEmpacar = c.hayCajasPorEmpacar;
            volEmpacado = 0.0;
            for (vector<PackedBox>::iterator p_it = empacados.begin(); p_it < empacados.end(); ++p_it) {
                volEmpacado += c.boxest[(*p_it).id].volumenReal;
            }
            volEmpacado /= (double)totalVolumen;
        }
    }

    // Métodos

private:
    void ActualizarMinD()
    {
        minDx = ContenedorDimx;
        minDy = ContenedorDimy;
        minDz = ContenedorDimz;
        for (vector<BoxType>::iterator b_it = boxest.begin(); b_it < boxest.end(); ++b_it) {
            if ((*b_it).cantidad > 0) {
                if ((*b_it).minDx < minDx) minDx = (*b_it).minDx;
                if ((*b_it).minDy < minDy) minDy = (*b_it).minDy;
                if ((*b_it).minDz < minDz) minDz = (*b_it).minDz;
            }
        }
    }

    // Espacios maximales

    void EliminarEspaciosMinD() {
        for (int i = (int)espacios.size() - 1; i >= 0; --i) {
            MaximalSpace& e = espacios[i];
            if (e.dx < minDx || e.dy < minDy || e.dz < minDz) espacios.erase(espacios.begin() + i);
        }
    }
    void DeterminarIndicesCajasEspacio(vector<MaximalSpace>::iterator const& e_it) {
        (*e_it).cambio = false;
        if ((*e_it).z1 > 0) {
            (*e_it).indCajas.clear();
            (*e_it).maxSoporte = 0;
            if (r_juntarEspacios == 0) {
                for (vector<PackedBox>::iterator p = empacados.begin(); p < empacados.end(); ++p) {
                    if ((*p).z2 == (*e_it).z1) {
                        if (!((*p).x1 >= (*e_it).x2 || (*e_it).x1 >= (*p).x2)) {
                            if (!((*p).y1 >= (*e_it).y2 || (*e_it).y1 >= (*p).y2)) {
                                (*e_it).indCajas.push_back(distance(empacados.begin(), p));
                                (*e_it).maxSoporte = max((*e_it).maxSoporte, (*p).presion);
                                return;
                            }
                        }
                    }
                }
            }
            else {
                for (vector<PackedBox>::iterator p = empacados.begin(); p < empacados.end(); ++p) {
                    if ((*p).z2 == (*e_it).z1) {
                        if (!((*p).x1 >= (*e_it).x2 || (*e_it).x1 >= (*p).x2)) {
                            if (!((*p).y1 >= (*e_it).y2 || (*e_it).y1 >= (*p).y2)) {
                                (*e_it).indCajas.push_back(distance(empacados.begin(), p));
                                (*e_it).maxSoporte = max((*e_it).maxSoporte, (*p).presion);
                            }
                        }
                    }
                }
            }
        }
    }
    void ActualizarIndicesCajasDeEspacios() {
        if (r_maxPresionItems && r_estabilidad == 1) { // Pueden haber espacios maximales que no estan soportados por ninguna caja y pueden ser eliminados
            for (int i = (int)espacios.size() - 1; i >= 0; --i) {
                vector<MaximalSpace>::iterator e_it = espacios.begin() + i;
                if ((*e_it).cambio) {
                    DeterminarIndicesCajasEspacio(e_it);
                    if ((*e_it).indCajas.size() == 0 && (*e_it).z1 > 0) {
                        espacios.erase(e_it);
                    }
                }
            }
        }
        else {
            for (vector<MaximalSpace>::iterator e_it = espacios.begin(); e_it < espacios.end(); ++e_it) {
                if ((*e_it).cambio) DeterminarIndicesCajasEspacio(e_it);
            }
        }
    }
    void JuntarEspacios(int desde) {
        //desde = 0;
        // Eliminar espacios contenidos

        for (int i = (int)espacios.size() - 1; i >= desde; --i) {
            MaximalSpace e2 = espacios[i];
            int viejosBorrados = 0;
            for (int j = i - 1; j >= 0; --j) {
                MaximalSpace e1 = espacios[j];

                // Ver si e1 contiene a e2

                if (e1.z1 <= e2.z1 && e2.z2 <= e1.z2) {
                    if (e1.y1 <= e2.y1 && e2.y2 <= e1.y2) {
                        if (e1.x1 <= e2.x1 && e2.x2 <= e1.x2) {
                            espacios.erase(espacios.begin() + (i - viejosBorrados));
                            break;
                        }
                    }
                }

                // Ver si e2 contiene a e1

                if (e2.z1 <= e1.z1 && e1.z2 <= e2.z2) {
                    if (e2.y1 <= e1.y1 && e1.y2 <= e2.y2) {
                        if (e2.x1 <= e1.x1 && e1.x2 <= e2.x2) {
                            espacios.erase(espacios.begin() + j);
                            ++viejosBorrados;
                            //break;
                        }
                    }
                }
            }
            i -= viejosBorrados;
        }

        // Juntar y expandir
        ++holi;
        if (r_juntarEspacios == 1) { // Se juntan espacios pero no se expanden
            if (r_estabilidad == 0) {
                bool seguir = true;
                while (seguir) {
                    seguir = false;
                    for (int i = (int)espacios.size() - 1; i >= desde; --i) {
                        MaximalSpace& e2 = espacios[i];
                        for (int j = i - 1; j >= 0; --j) {
                            MaximalSpace& e1 = espacios[j];

                            // Ver si se pueden juntar

                            if (e1.z1 == e2.z1 && e1.z2 == e2.z2) {
                                if (e1.y1 == e2.y1 && e1.y2 == e2.y2) {
                                    if (!(e1.x1 > e2.x2 || e2.x1 > e1.x2)) {
                                        e1.x1 = min(e1.x1, e2.x1);
                                        e1.x2 = max(e1.x2, e2.x2);
                                        e1.dx = e1.x2 - e1.x1;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }
                                else if (e1.x1 == e2.x1 && e1.x2 == e2.x2) {
                                    if (!(e1.y1 > e2.y2 || e2.y1 > e1.y2)) {
                                        e1.y1 = min(e1.y1, e2.y1);
                                        e1.y2 = max(e1.y2, e2.y2);
                                        e1.dy = e1.y2 - e1.y1;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }
                            }
                            else if (e1.x1 == e2.x1 && e1.x2 == e2.x2) {
                                if (e1.y1 == e2.y1 && e1.y2 == e2.y2) {
                                    if (!(e1.z1 > e2.z2 || e2.z1 > e1.z2)) {
                                        e1.z1 = min(e1.z1, e2.z1);
                                        e1.z2 = max(e1.z2, e2.z2);
                                        e1.dz = e1.z2 - e1.z1;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else {
                bool seguir = true;
                while (seguir) {
                    seguir = false;
                    for (int i = (int)espacios.size() - 1; i >= desde; --i) {
                        MaximalSpace& e2 = espacios[i];
                        for (int j = i - 1; j >= 0; --j) {
                            MaximalSpace& e1 = espacios[j];

                            // Ver si se pueden juntar

                            if (e1.z1 == e2.z1) {
                                if (e1.y1 == e2.y1 && e1.y2 == e2.y2) {
                                    if (!(e1.x1 > e2.x2 || e2.x1 > e1.x2)) {
                                        e1.x1 = min(e1.x1, e2.x1);
                                        e1.x2 = max(e1.x2, e2.x2);
                                        e1.dx = e1.x2 - e1.x1;
                                        e1.cambio = true;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }
                                else if (e1.x1 == e2.x1 && e1.x2 == e2.x2) {
                                    if (!(e1.y1 > e2.y2 || e2.y1 > e1.y2)) {
                                        e1.y1 = min(e1.y1, e2.y1);
                                        e1.y2 = max(e1.y2, e2.y2);
                                        e1.dy = e1.y2 - e1.y1;
                                        e1.cambio = true;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        else if (r_juntarEspacios == 2) { // Se juntan y se expanden espacios
            if (r_estabilidad == 0) {
                bool seguir = true;
                while (seguir) {
                    seguir = false;
                    for (int i = (int)espacios.size() - 1; i >= desde; --i) {
                        MaximalSpace& e2 = espacios[i];
                        for (int j = i - 1; j >= 0; --j) {
                            MaximalSpace& e1 = espacios[j];

                            // Ver si se pueden juntar

                            if (e1.z1 == e2.z1 && e1.z2 == e2.z2) {
                                if (e1.y1 == e2.y1 && e1.y2 == e2.y2) {
                                    if (!(e1.x1 > e2.x2 || e2.x1 > e1.x2)) {
                                        e1.x1 = min(e1.x1, e2.x1);
                                        e1.x2 = max(e1.x2, e2.x2);
                                        e1.dx = e1.x2 - e1.x1;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }
                                else if (e1.x1 == e2.x1 && e1.x2 == e2.x2) {
                                    if (!(e1.y1 > e2.y2 || e2.y1 > e1.y2)) {
                                        e1.y1 = min(e1.y1, e2.y1);
                                        e1.y2 = max(e1.y2, e2.y2);
                                        e1.dy = e1.y2 - e1.y1;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }
                            }
                            else if (e1.x1 == e2.x1 && e1.x2 == e2.x2) {
                                if (e1.y1 == e2.y1 && e1.y2 == e2.y2) {
                                    if (!(e1.z1 > e2.z2 || e2.z1 > e1.z2)) {
                                        e1.z1 = min(e1.z1, e2.z1);
                                        e1.z2 = max(e1.z2, e2.z2);
                                        e1.dz = e1.z2 - e1.z1;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }
                            }

                            // Ver si e1 se puede expandir con e2

                            if (e2.z1 <= e1.z1 && e1.z2 <= e2.z2) {
                                if (e2.y1 <= e1.y1 && e1.y2 <= e2.y2) {
                                    if (!(e1.x1 > e2.x2 || e2.x1 > e1.x2)) {
                                        int mm = min(e1.x1, e2.x1);
                                        int MM = max(e1.x2, e2.x2);
                                        if (e1.x1 != mm || e1.x2 != MM) {
                                            e1.x1 = mm;
                                            e1.x2 = MM;
                                            e1.dx = MM - mm;
                                            seguir = true;
                                        }
                                    }
                                }
                                if (e2.x1 <= e1.x1 && e1.x2 <= e2.x2) {
                                    if (!(e1.y1 > e2.y2 || e2.y1 > e1.y2)) {
                                        int mm = min(e1.y1, e2.y1);
                                        int MM = max(e1.y2, e2.y2);
                                        if (e1.y1 != mm || e1.y2 != MM) {
                                            e1.y1 = mm;
                                            e1.y2 = MM;
                                            e1.dy = MM - mm;
                                            seguir = true;
                                        }
                                    }
                                }
                            }
                            if (e2.y1 <= e1.y1 && e1.y2 <= e2.y2) {
                                if (e2.x1 <= e1.x1 && e1.x2 <= e2.x2) {
                                    if (!(e1.z1 > e2.z2 || e2.z1 > e1.z2)) {
                                        int mm = min(e1.z1, e2.z1);
                                        int MM = max(e1.z2, e2.z2);
                                        if (e1.z1 != mm || e1.z2 != MM) {
                                            e1.z1 = mm;
                                            e1.z2 = MM;
                                            e1.dz = MM - mm;
                                            seguir = true;
                                        }
                                    }
                                }
                            }

                            // Ver si e2 se puede expandir con e1

                            if (e1.z1 <= e2.z1 && e2.z2 <= e1.z2) {
                                if (e1.y1 <= e2.y1 && e2.y2 <= e1.y2) {
                                    if (!(e2.x1 > e1.x2 || e1.x1 > e2.x2)) {
                                        int mm = min(e2.x1, e1.x1);
                                        int MM = max(e2.x2, e1.x2);
                                        if (e2.x1 != mm || e2.x2 != MM) {
                                            e2.x1 = mm;
                                            e2.x2 = MM;
                                            e2.dx = MM - mm;
                                            seguir = true;
                                        }
                                    }
                                }
                                if (e1.x1 <= e2.x1 && e2.x2 <= e1.x2) {
                                    if (!(e2.y1 > e1.y2 || e1.y1 > e2.y2)) {
                                        int mm = min(e2.y1, e1.y1);
                                        int MM = max(e2.y2, e1.y2);
                                        if (e2.y1 != mm || e2.y2 != MM) {
                                            e2.y1 = mm;
                                            e2.y2 = MM;
                                            e2.dy = MM - mm;
                                            seguir = true;
                                        }
                                    }
                                }
                            }
                            if (e1.y1 <= e2.y1 && e2.y2 <= e1.y2) {
                                if (e1.x1 <= e2.x1 && e2.x2 <= e1.x2) {
                                    if (!(e2.z1 > e1.z2 || e1.z1 > e2.z2)) {
                                        int mm = min(e2.z1, e1.z1);
                                        int MM = max(e2.z2, e1.z2);
                                        if (e2.z1 != mm || e2.z2 != MM) {
                                            e2.z1 = mm;
                                            e2.z2 = MM;
                                            e2.dz = MM - mm;
                                            seguir = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else if (r_estabilidad == 1) {
                bool seguir = true;
                while (seguir) {
                    seguir = false;
                    for (int i = (int)espacios.size() - 1; i >= desde; --i) {
                        MaximalSpace& e2 = espacios[i];
                        for (int j = i - 1; j >= 0; --j) {
                            MaximalSpace& e1 = espacios[j];

                            // Ver si se pueden juntar

                            if (e1.z1 == e2.z1) {
                                if (e1.y1 == e2.y1 && e1.y2 == e2.y2) {
                                    if (!(e1.x1 > e2.x2 || e2.x1 > e1.x2)) {
                                        e1.x1 = min(e1.x1, e2.x1);
                                        e1.x2 = max(e1.x2, e2.x2);
                                        e1.dx = e1.x2 - e1.x1;
                                        e1.cambio = true;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }
                                else if (e1.x1 == e2.x1 && e1.x2 == e2.x2) {
                                    if (!(e1.y1 > e2.y2 || e2.y1 > e1.y2)) {
                                        e1.y1 = min(e1.y1, e2.y1);
                                        e1.y2 = max(e1.y2, e2.y2);
                                        e1.dy = e1.y2 - e1.y1;
                                        e1.cambio = true;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }
                            }

                            // Ver si e1 se puede expandir con e2

                            if (e2.z1 <= e1.z1) {
                                if (e2.y1 <= e1.y1 && e1.y2 <= e2.y2) {
                                    if (!(e1.x1 > e2.x2 || e2.x1 > e1.x2)) {
                                        int mm = min(e1.x1, e2.x1);
                                        int MM = max(e1.x2, e2.x2);
                                        if (e1.x1 != mm || e1.x2 != MM) {
                                            e1.x1 = mm;
                                            e1.x2 = MM;
                                            e1.dx = MM - mm;
                                            e1.cambio = true;
                                            seguir = true;
                                        }
                                    }
                                }
                                if (e2.x1 <= e1.x1 && e1.x2 <= e2.x2) {
                                    if (!(e1.y1 > e2.y2 || e2.y1 > e1.y2)) {
                                        int mm = min(e1.y1, e2.y1);
                                        int MM = max(e1.y2, e2.y2);
                                        if (e1.y1 != mm || e1.y2 != MM) {
                                            e1.y1 = mm;
                                            e1.y2 = MM;
                                            e1.dy = MM - mm;
                                            e1.cambio = true;
                                            seguir = true;
                                        }
                                    }
                                }
                            }

                            // Ver si e2 se puede expandir con e1

                            if (e1.z1 <= e2.z1) {
                                if (e1.y1 <= e2.y1 && e2.y2 <= e1.y2) {
                                    if (!(e2.x1 > e1.x2 || e1.x1 > e2.x2)) {
                                        int mm = min(e2.x1, e1.x1);
                                        int MM = max(e2.x2, e1.x2);
                                        if (e2.x1 != mm || e2.x2 != MM) {
                                            e2.x1 = mm;
                                            e2.x2 = MM;
                                            e2.dx = MM - mm;
                                            e2.cambio = true;
                                            seguir = true;
                                        }
                                    }
                                }
                                if (e1.x1 <= e2.x1 && e2.x2 <= e1.x2) {
                                    if (!(e2.y1 > e1.y2 || e1.y1 > e2.y2)) {
                                        int mm = min(e2.y1, e1.y1);
                                        int MM = max(e2.y2, e1.y2);
                                        if (e2.y1 != mm || e2.y2 != MM) {
                                            e2.y1 = mm;
                                            e2.y2 = MM;
                                            e2.dy = MM - mm;
                                            e2.cambio = true;
                                            seguir = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else if (r_estabilidad == 2) {
                bool seguir = true;
                while (seguir) {
                    seguir = false;
                    for (int i = (int)espacios.size() - 1; i >= desde; --i) {
                        MaximalSpace& e2 = espacios[i];
                        for (int j = i - 1; j >= 0; --j) {
                            MaximalSpace& e1 = espacios[j];

                            // Ver si se pueden juntar

                            if (e1.z1 == e2.z1 && e1.z2 == e2.z2) {
                                if (e1.y1 == e2.y1 && e1.y2 == e2.y2) {
                                    if (!(e1.x1 > e2.x2 || e2.x1 > e1.x2)) {
                                        e1.x1 = min(e1.x1, e2.x1);
                                        e1.x2 = max(e1.x2, e2.x2);
                                        e1.dx = e1.x2 - e1.x1;
                                        e1.cambio = true;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }
                                else if (e1.x1 == e2.x1 && e1.x2 == e2.x2) {
                                    if (!(e1.y1 > e2.y2 || e2.y1 > e1.y2)) {
                                        e1.y1 = min(e1.y1, e2.y1);
                                        e1.y2 = max(e1.y2, e2.y2);
                                        e1.dy = e1.y2 - e1.y1;
                                        e1.cambio = true;
                                        espacios.erase(espacios.begin() + i);
                                        seguir = true;
                                        break;
                                    }
                                }

                                // Ver si e1 se puede expandir con e2

                                if (e2.y1 <= e1.y1 && e1.y2 <= e2.y2) {
                                    if (!(e1.x1 > e2.x2 || e2.x1 > e1.x2)) {
                                        int mm = min(e1.x1, e2.x1);
                                        int MM = max(e1.x2, e2.x2);
                                        if (e1.x1 != mm || e1.x2 != MM) {
                                            e1.x1 = mm;
                                            e1.x2 = MM;
                                            e1.dx = MM - mm;
                                            e1.cambio = true;
                                            seguir = true;
                                        }
                                    }
                                }
                                if (e2.x1 <= e1.x1 && e1.x2 <= e2.x2) {
                                    if (!(e1.y1 > e2.y2 || e2.y1 > e1.y2)) {
                                        int mm = min(e1.y1, e2.y1);
                                        int MM = max(e1.y2, e2.y2);
                                        if (e1.y1 != mm || e1.y2 != MM) {
                                            e1.y1 = mm;
                                            e1.y2 = MM;
                                            e1.dy = MM - mm;
                                            e1.cambio = true;
                                            seguir = true;
                                        }
                                    }
                                }

                                // Ver si e2 se puede expandir con e1

                                if (e1.y1 <= e2.y1 && e2.y2 <= e1.y2) {
                                    if (!(e2.x1 > e1.x2 || e1.x1 > e2.x2)) {
                                        int mm = min(e2.x1, e1.x1);
                                        int MM = max(e2.x2, e1.x2);
                                        if (e2.x1 != mm || e2.x2 != MM) {
                                            e2.x1 = mm;
                                            e2.x2 = MM;
                                            e2.dx = MM - mm;
                                            e2.cambio = true;
                                            seguir = true;
                                        }
                                    }
                                }
                                if (e1.x1 <= e2.x1 && e2.x2 <= e1.x2) {
                                    if (!(e2.y1 > e1.y2 || e1.y1 > e2.y2)) {
                                        int mm = min(e2.y1, e1.y1);
                                        int MM = max(e2.y2, e1.y2);
                                        if (e2.y1 != mm || e2.y2 != MM) {
                                            e2.y1 = mm;
                                            e2.y2 = MM;
                                            e2.dy = MM - mm;
                                            e2.cambio = true;
                                            seguir = true;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    int ActualizarEspacios(Bloque const& c) {
        int resp = espacios.size();
        if (r_estabilidad == 0) {
            for (int i = resp - 1; i >= 0; --i) {
                MaximalSpace e = espacios[i];
                if (!(e.x1 >= c.x2 || e.y1 >= c.y2 || e.z1 >= c.z2 || c.x1 >= e.x2 || c.y1 >= e.y2 || c.z1 >= e.z2)) {

                    // Espacios en x

                    if (c.x1 > e.x1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, c.x1, e.y2, e.z2));
                    if (e.x2 > c.x2) espacios.push_back(MaximalSpace(c.x2, e.y1, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en y

                    if (c.y1 > e.y1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, e.x2, c.y1, e.z2));
                    if (e.y2 > c.y2) espacios.push_back(MaximalSpace(e.x1, c.y2, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en z

                    if (c.z1 > e.z1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, e.x2, e.y2, c.z1));
                    if (e.z2 > c.z2) espacios.push_back(MaximalSpace(e.x1, e.y1, c.z2, e.x2, e.y2, e.z2));
                    espacios.erase(espacios.begin() + i);
                    --resp;
                }
            }
        }
        else if (r_estabilidad == 1) {
            for (int i = resp - 1; i >= 0; --i) {
                MaximalSpace e = espacios[i];
                if (!(e.x1 >= c.x2 || e.y1 >= c.y2 || e.z1 >= c.z2 || c.x1 >= e.x2 || c.y1 >= e.y2)) {

                    // Espacios en x

                    if (c.x1 > e.x1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, c.x1, e.y2, e.z2));
                    if (e.x2 > c.x2) espacios.push_back(MaximalSpace(c.x2, e.y1, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en y

                    if (c.y1 > e.y1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, e.x2, c.y1, e.z2));
                    if (e.y2 > c.y2) espacios.push_back(MaximalSpace(e.x1, c.y2, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en z

                    if (e.z2 > c.z2) espacios.push_back(MaximalSpace(e.x1, e.y1, c.z2, e.x2, e.y2, e.z2));
                    espacios.erase(espacios.begin() + i);
                    --resp;
                }
            }
        }
        else {
            for (int i = resp - 1; i >= 0; --i) {
                MaximalSpace e = espacios[i];
                if (!(e.x1 >= c.x2 || e.y1 >= c.y2 || e.z1 >= c.z2 || c.x1 >= e.x2 || c.y1 >= e.y2)) {

                    // Espacios en x

                    if (c.x1 > e.x1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, c.x1, e.y2, e.z2));
                    if (e.x2 > c.x2) espacios.push_back(MaximalSpace(c.x2, e.y1, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en y

                    if (c.y1 > e.y1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, e.x2, c.y1, e.z2));
                    if (e.y2 > c.y2) espacios.push_back(MaximalSpace(e.x1, c.y2, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en z

                    if (e.z2 > c.z2) espacios.push_back(MaximalSpace(c.x1, c.y1, c.z2, c.x2, c.y2, e.z2));
                    espacios.erase(espacios.begin() + i);
                    --resp;
                }
            }
        }
        return resp;
    }
    int ActualizarEspacios(PackedBox const& c) {
        int resp = espacios.size();
        if (r_estabilidad == 0) {
            for (int i = resp - 1; i >= 0; --i) {
                MaximalSpace e = espacios[i];
                if (!(e.x1 >= c.x2 || e.y1 >= c.y2 || e.z1 >= c.z2 || c.x1 >= e.x2 || c.y1 >= e.y2 || c.z1 >= e.z2)) {

                    // Espacios en x

                    if (c.x1 > e.x1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, c.x1, e.y2, e.z2));
                    if (e.x2 > c.x2) espacios.push_back(MaximalSpace(c.x2, e.y1, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en y

                    if (c.y1 > e.y1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, e.x2, c.y1, e.z2));
                    if (e.y2 > c.y2) espacios.push_back(MaximalSpace(e.x1, c.y2, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en z

                    if (c.z1 > e.z1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, e.x2, e.y2, c.z1));
                    if (e.z2 > c.z2) espacios.push_back(MaximalSpace(e.x1, e.y1, c.z2, e.x2, e.y2, e.z2));
                    espacios.erase(espacios.begin() + i);
                    --resp;
                }
            }
        }
        else if (r_estabilidad == 1) {
            for (int i = resp - 1; i >= 0; --i) {
                MaximalSpace e = espacios[i];
                if (!(e.x1 >= c.x2 || e.y1 >= c.y2 || e.z1 >= c.z2 || c.x1 >= e.x2 || c.y1 >= e.y2)) {

                    // Espacios en x

                    if (c.x1 > e.x1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, c.x1, e.y2, e.z2));
                    if (e.x2 > c.x2) espacios.push_back(MaximalSpace(c.x2, e.y1, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en y

                    if (c.y1 > e.y1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, e.x2, c.y1, e.z2));
                    if (e.y2 > c.y2) espacios.push_back(MaximalSpace(e.x1, c.y2, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en z

                    if (e.z2 > c.z2) espacios.push_back(MaximalSpace(e.x1, e.y1, c.z2, e.x2, e.y2, e.z2));
                    espacios.erase(espacios.begin() + i);
                    --resp;
                }
            }
        }
        else {
            for (int i = resp - 1; i >= 0; --i) {
                MaximalSpace e = espacios[i];
                if (!(e.x1 >= c.x2 || e.y1 >= c.y2 || e.z1 >= c.z2 || c.x1 >= e.x2 || c.y1 >= e.y2)) {

                    // Espacios en x

                    if (c.x1 > e.x1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, c.x1, e.y2, e.z2));
                    if (e.x2 > c.x2) espacios.push_back(MaximalSpace(c.x2, e.y1, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en y

                    if (c.y1 > e.y1) espacios.push_back(MaximalSpace(e.x1, e.y1, e.z1, e.x2, c.y1, e.z2));
                    if (e.y2 > c.y2) espacios.push_back(MaximalSpace(e.x1, c.y2, e.z1, e.x2, e.y2, e.z2));

                    // Espacios en z

                    if (e.z2 > c.z2) espacios.push_back(MaximalSpace(c.x1, c.y1, c.z2, c.x2, c.y2, e.z2));
                    espacios.erase(espacios.begin() + i);
                    --resp;
                }
            }
        }
        return resp;
    }
    void ActualizarEspaciosCapas(int const& x0) {
        if (r_estabilidad == 0) {
            for (int i = (int)espacios.size() - 1; i >= 0; --i) {
                MaximalSpace& e = espacios[i];
                if (e.x2 <= x0) espacios.erase(espacios.begin() + i);
                else {
                    e.x1 = x0;
                    e.dx = e.x2 - e.x1;
                    e.DeterminarEsquina();
                }
            }
        }
        else {
            for (int i = (int)espacios.size() - 1; i >= 0; --i) {
                MaximalSpace& e = espacios[i];
                if (e.x2 <= x0) espacios.erase(espacios.begin() + i);
                else {
                    e.x1 = x0;
                    e.dx = e.x2 - e.x1;
                    e.cambio = true;
                    e.DeterminarEsquina();
                }
            }
        }
    }
    void ActualizarEspaciosVisibilidad() {
        if (r_estabilidad == 0) {
            for (int i = (int)espacios.size() - 1; i >= 0; --i) {
                if (espacios[i].x2 < ContenedorDimx) {
                    espacios.erase(espacios.begin() + i);
                }
            }
        }
        else {
            for (int i = (int)espacios.size() - 1; i >= 0; --i) {
                MaximalSpace e = espacios[i];
                if (e.x2 < ContenedorDimx) {
                    for (vector<PackedBox>::reverse_iterator p_it = empacados.rbegin(); p_it != empacados.rend(); ++p_it) {
                        if (e.x1 < (*p_it).x2) {
                            if (!(e.z1 >= (*p_it).z2 || (*p_it).z1 >= e.z2)) {
                                if (!(e.y1 >= (*p_it).y2 || (*p_it).y1 >= e.y2)) {
                                    espacios.erase(espacios.begin() + i);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // Selección Espacios Maximales

    void SeleccionarEspacio_Fondo_Altura(vector<MaximalSpace>::iterator& MejorEspacio) {
        //sort(espacios.begin(), espacios.end(), Orden_Fondo_Altura);
        //uniform_int_distribution<int> distri(0, miBusquedaBinaria(espacios.front().x1 + (double)(espacios.back().x1 - espacios.front().x1) * alphaEM + error));
        //MejorEspacio = espacios.begin() + distri(generator);
        //MejorEspacio = espacios.begin();
        for (vector<MaximalSpace>::iterator e_it = espacios.begin() + 1; e_it < espacios.end(); ++e_it) {
            if ((*e_it).x1 < (*MejorEspacio).x1) MejorEspacio = e_it;
            else if ((*e_it).x1 == (*MejorEspacio).x1) {
                if ((*e_it).z1 > (*MejorEspacio).z1) MejorEspacio = e_it;
                else if ((*e_it).z1 == (*MejorEspacio).z1) {
                    if ((*e_it).metrica < (*MejorEspacio).metrica) MejorEspacio = e_it;
                }
            }
        }
    }
    void SeleccionarEspacio_Fondo_Bajo(vector<MaximalSpace>::iterator& MejorEspacio) {
        //sort(espacios.begin(), espacios.end(), Orden_Fondo_Bajo);
        //uniform_int_distribution<int> distri(0, miBusquedaBinaria(espacios.front().x1 + (double)(espacios.back().x1 - espacios.front().x1) * alphaEM + error));
        //MejorEspacio = espacios.begin() + distri(generator);
        //MejorEspacio = espacios.begin();
        for (vector<MaximalSpace>::iterator e_it = espacios.begin() + 1; e_it < espacios.end(); ++e_it) {
            if ((*e_it).x1 < (*MejorEspacio).x1) MejorEspacio = e_it;
            else if ((*e_it).x1 == (*MejorEspacio).x1) {
                if ((*e_it).z1 < (*MejorEspacio).z1) MejorEspacio = e_it;
                else if ((*e_it).z1 == (*MejorEspacio).z1) {
                    if ((*e_it).metrica < (*MejorEspacio).metrica) MejorEspacio = e_it;
                }
            }
        }
    }
    void SeleccionarEspacio_Altura_Fondo(vector<MaximalSpace>::iterator& MejorEspacio) {
        //sort(espacios.begin(), espacios.end(), Orden_Altura_Fondo);
        //uniform_int_distribution<int> distri(0, miBusquedaBinaria(espacios.front().z1 + (double)(espacios.back().z1 - espacios.front().z1) * alphaEM + error));
        //MejorEspacio = espacios.begin() + distri(generator);
        //MejorEspacio = espacios.begin();
        for (vector<MaximalSpace>::iterator e_it = espacios.begin() + 1; e_it < espacios.end(); ++e_it) {
            if ((*e_it).z1 > (*MejorEspacio).z1) MejorEspacio = e_it;
            else if ((*e_it).z1 == (*MejorEspacio).z1) {
                if ((*e_it).x1 < (*MejorEspacio).x1) MejorEspacio = e_it;
                else if ((*e_it).x1 == (*MejorEspacio).x1) {
                    if ((*e_it).metrica < (*MejorEspacio).metrica) MejorEspacio = e_it;
                }
            }
        }
    }

    // Selección Items

    int miBusquedaBinariaBF(vector<Bloque> const& b, double const& val) {
        if (b.back().bf1 < val) return b.size() - 1;
        int l = 0;
        int r = b.size() - 1;
        while (l < r) {
            int m = (l + r) / 2;
            if (b[m].bf1 > val) r = m - 1;
            else l = m + 1;
        }
        return r;
    }
    int miBusquedaBinariaVol(vector<Bloque> const& b, double const& val) {
        if (b.back().vol < val) return b.size() - 1;
        int l = 0;
        int r = b.size() - 1;
        while (l < r) {
            int m = (l + r) / 2;
            if (b[m].vol > val) r = m - 1;
            else l = m + 1;
        }
        return r;
    }
    bool BloqueSoportadoPoligono(Bloque& b) {
        if (b.z1 > 0) {
            for (vector<BloqueAux>::iterator b_it = b.cajas.begin(); b_it < b.cajas.end(); ++b_it) {
                if ((*b_it).z1 == b.z1) { // solo se verifica las cajas de abajo
                    double centroMasax = (double)((*b_it).x1 + (double)b.b.x / 2.0);
                    double centroMasay = (double)((*b_it).y1 + (double)b.b.y / 2.0);
                    vector<Point> intersecciones;
                    intersecciones.reserve(4 * b.indPB_EM.size());
                    for (vector<int>::iterator idPB = b.indPB_EM.begin(); idPB < b.indPB_EM.end(); ++idPB) {
                        PackedBox& p = empacados[*idPB];
                        if (!((*b_it).x1 >= p.x2 || p.x1 >= (*b_it).x1 + b.b.x)) {
                            if (!((*b_it).y1 >= p.y2 || p.y1 >= (*b_it).y1 + b.b.y)) {

                                // Se encuentran los puntos de intersección

                                vector<int> misx({ (*b_it).x1, (*b_it).x1 + b.b.x, p.x1, p.x2 });
                                sort(misx.begin(), misx.end());
                                vector<int> misy({ (*b_it).y1, (*b_it).y1 + b.b.y, p.y1, p.y2 });
                                sort(misy.begin(), misy.end());
                                intersecciones.push_back(Point(misx[1], misy[1]));
                                intersecciones.push_back(Point(misx[2], misy[1]));
                                intersecciones.push_back(Point(misx[1], misy[2]));
                                intersecciones.push_back(Point(misx[2], misy[2]));
                            }
                        }
                    }

                    // Determinar polígono

                    if (intersecciones.size() > 0) {
                        Poligono poligono(intersecciones);
                        if (!poligono.puntoDentroPoligono(centroMasax, centroMasay)) {
                            return false;
                        }
                    }
                    else return false;
                }
            }
        }
        return true;
    }
    bool BloqueSoportadoPresion(Bloque& b) {

        // Generar las piezas empacadas y se inicializan los soportes de las cajas más arriba

        vector<PackedBox> copias;
        copias.reserve(empacados.size() + b.cajas.size());
        for (vector<BloqueAux>::iterator b_it = b.cajas.begin(); b_it < b.cajas.end(); ++b_it) {
            if ((*b_it).z1 == b.z1) {
                if (b.z1 == 0) {
                    copias.push_back(PackedBox(copias.size() + empacados.size(), b, *b_it));
                }
                else {
                    double centroMasax = (double)((*b_it).x1 + (double)b.b.x / 2.0);
                    double centroMasay = (double)((*b_it).y1 + (double)b.b.y / 2.0);
                    vector<Point> intersecciones;
                    intersecciones.reserve(4 * b.indPB_EM.size());
                    vector<Soporte> misSoportes; misSoportes.reserve(b.indPB_EM.size());
                    for (vector<int>::iterator idPB = b.indPB_EM.begin(); idPB < b.indPB_EM.end(); ++idPB) {
                        PackedBox& p = empacados[*idPB];
                        if (!((*b_it).x1 >= p.x2 || p.x1 >= (*b_it).x1 + b.b.x)) {
                            if (!((*b_it).y1 >= p.y2 || p.y1 >= (*b_it).y1 + b.b.y)) {

                                // Se encuentran los puntos de intersección

                                vector<int> misx({ (*b_it).x1, (*b_it).x1 + b.b.x, p.x1, p.x2 });
                                sort(misx.begin(), misx.end());
                                vector<int> misy({ (*b_it).y1, (*b_it).y1 + b.b.y, p.y1, p.y2 });
                                sort(misy.begin(), misy.end());
                                intersecciones.push_back(Point(misx[1], misy[1]));
                                intersecciones.push_back(Point(misx[2], misy[1]));
                                intersecciones.push_back(Point(misx[1], misy[2]));
                                intersecciones.push_back(Point(misx[2], misy[2]));

                                // Se crea el soporte

                                misSoportes.push_back(Soporte(*idPB, (misx[2] - misx[1]) * (misy[2] - misy[1]), (double)(misx[2] + misx[1]) / 2.0f, (double)(misy[2] + misy[1]) / 2.0f, p.presion));
                            }
                        }
                    }

                    // Determinar polígono

                    if (intersecciones.size() > 0) {
                        Poligono poligono(intersecciones);
                        if (!poligono.puntoDentroPoligono(centroMasax, centroMasay)) {
                            return false;
                        }
                        copias.push_back(PackedBox(copias.size() + empacados.size(), b, *b_it, centroMasax, centroMasay, poligono, misSoportes));
                        if (copias.back().z2 == b.z2 && copias.back().z1 > 0) {
                            if (!copias.back().ActualizarSoportes()) {
                                return false;
                            }
                        }
                    }
                    else return false;
                }
            }
            else {
                double centroMasax = (double)((*b_it).x1 + (double)b.b.x / 2.0);
                double centroMasay = (double)((*b_it).y1 + (double)b.b.y / 2.0);
                vector<Point> intersecciones;
                intersecciones.reserve(4 * b.indPB_EM.size());
                vector<Soporte> misSoportes; misSoportes.reserve(b.indPB_EM.size());
                for (vector<PackedBox>::iterator p_it = copias.begin(); p_it < copias.end(); ++p_it) {
                    if ((*p_it).z2 == (*b_it).z1) {
                        if (!((*b_it).x1 >= (*p_it).x2 || (*p_it).x1 >= (*b_it).x1 + b.b.x)) {
                            if (!((*b_it).y1 >= (*p_it).y2 || (*p_it).y1 >= (*b_it).y1 + b.b.y)) {

                                // Se encuentran los puntos de intersección

                                vector<int> misx({ (*b_it).x1, (*b_it).x1 + b.b.x, (*p_it).x1, (*p_it).x2 });
                                sort(misx.begin(), misx.end());
                                vector<int> misy({ (*b_it).y1, (*b_it).y1 + b.b.y, (*p_it).y1, (*p_it).y2 });
                                sort(misy.begin(), misy.end());
                                intersecciones.push_back(Point(misx[1], misy[1]));
                                intersecciones.push_back(Point(misx[2], misy[1]));
                                intersecciones.push_back(Point(misx[1], misy[2]));
                                intersecciones.push_back(Point(misx[2], misy[2]));

                                // Se crea el soporte

                                misSoportes.push_back(Soporte((*p_it).ind, (misx[2] - misx[1]) * (misy[2] - misy[1]), (double)(misx[2] + misx[1]) / 2.0f, (double)(misy[2] + misy[1]) / 2.0f, (*p_it).presion));
                            }
                        }
                    }
                }

                // Determinar polígono

                if (intersecciones.size() > 0) {
                    Poligono poligono(intersecciones);
                    copias.push_back(PackedBox(copias.size(), b, *b_it, centroMasax, centroMasay, poligono, misSoportes));
                    if (copias.back().z2 == b.z2 && copias.back().z1 > 0) {
                        if (!copias.back().ActualizarSoportes()) {
                            return false;
                        }
                    }
                }
                else return false;
            }
        }

        // Determinar arbol de actualización

        if (b.z1 > 0) {
            int hasta = copias.size();
            for (int i = 0; i < hasta; ++i) {
                PackedBox c = copias[i];
                if (c.z1 == b.z1) {
                    int siguiente = copias.size();
                    while (c.z1 > 0) {
                        for (vector<Soporte>::iterator s_it = c.soportes.begin(); s_it < c.soportes.end(); ++s_it) {
                            bool nuevaPieza = true;
                            for (vector<PackedBox>::iterator p_it = copias.begin() + hasta; p_it < copias.end(); ++p_it) {
                                if ((*s_it).indPB == (*p_it).ind) {
                                    nuevaPieza = false;
                                    break;
                                }
                            }
                            if (nuevaPieza) {
                                copias.push_back(empacados[(*s_it).indPB]);
                            }
                        }
                        if (siguiente < copias.size()) {
                            for (; siguiente < copias.size(); ++siguiente) {
                                c = copias[siguiente];
                                if (c.z1 > 0) break;
                            }
                            ++siguiente;
                        }
                        else break;
                    }
                }
            }
        }
        sort(copias.begin(), copias.end(), [](PackedBox const& p1, PackedBox const& p2)->bool {return p1.z1 > p2.z1; });

        // Actualizar árbol

        int z1Ini = copias.front().z1;
        for (vector<PackedBox>::iterator p_it = copias.begin(); p_it < copias.end(); ++p_it) {
            if (z1Ini == (*p_it).z1) {
                for (vector<Soporte>::iterator s_it = (*p_it).soportes.begin(); s_it < (*p_it).soportes.end(); ++s_it) {

                    // Actualizar soporte

                    vector<PackedBox>::iterator p_it2 = p_it + 1;
                    for (; p_it2 < copias.end(); ++p_it2) {
                        if ((*p_it2).ind == (*s_it).indPB) {
                            (*p_it2).ActualizarCentroDeMasa(*s_it);
                            break;
                        }
                    }
                }
            }
            else {
                for (vector<PackedBox>::iterator p_it2 = p_it; p_it2 < copias.end(); ++p_it2) {
                    if ((*p_it2).z2 == z1Ini && (*p_it2).z1 > 0) {
                        if ((*p_it2).centroMasaEnPoligono()) {
                            if (!(*p_it2).ActualizarSoportes()) {
                                return false;
                            }
                        }
                        else {
                            return false;
                        }
                    }
                }
                z1Ini = (*p_it).z1;
                --p_it;
            }
        }
        return true;
    }
    bool EsBloqueNuevo(vector<Bloque>& opciones, int const& desde, Bloque nuevo) {
        for (vector<Bloque>::iterator b_it = opciones.begin() + desde; b_it < opciones.end(); ++b_it) {
            if ((*b_it) == nuevo) return false;
        }
        return true;
    }
    void MeterBloqueNuevo(vector<Bloque>& opciones, int const& desde, Bloque nuevo) {
        for (vector<Bloque>::iterator b_it = opciones.begin() + desde; b_it < opciones.end(); ++b_it) {
            if ((*b_it) == nuevo) return;
        }
        opciones.push_back(nuevo);
    }
    void GenerarOpciones(vector<Bloque>& opciones, MaximalSpace const& e) {
        opciones.reserve(maxNBloques);
        if (r_estabilidad == 1) {
            if (dist01(generator) < 0.5) { // Bloques de 1 sola caja en x
                for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                    if ((*bt).cantidad == 1) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                Bloque nuevoB(*bt, *b, e);
                                if (BloqueSoportadoPoligono(nuevoB)) {
                                    opciones.push_back(nuevoB);
                                }
                            }
                        }
                    }
                    else if ((*bt).cantidad > 1) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                int desde = opciones.size();

                                // y z

                                int ny = min(e.dy / (*b).y, (*bt).cantidad);
                                while (ny > 0) {
                                    Bloque nuevo(*bt, *b, e, 1, ny, 1);
                                    if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                        if (BloqueSoportadoPoligono(nuevo)) {
                                            opciones.push_back(nuevo);
                                            int nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                            while (nz > 1) {
                                                nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        break;
                                                    }
                                                }
                                                --nz;
                                            }
                                            break;
                                        }
                                    }
                                    --ny;
                                }

                                // z y

                                int nz = min(e.dz / (*b).z, (*bt).cantidad);
                                while (nz > 0) {
                                    Bloque nuevo(*bt, *b, e, 1, 1, nz);
                                    if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                        if (BloqueSoportadoPoligono(nuevo)) {
                                            opciones.push_back(nuevo);
                                            int ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                            while (ny > 1) {
                                                nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        break;
                                                    }
                                                }
                                                --ny;
                                            }
                                            break;
                                        }
                                    }
                                    --nz;
                                }
                            }
                        }
                    }
                }
            }
            else { // Bloques de varias cajas en x
                for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                    if ((*bt).cantidad == 1) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                Bloque nuevoB(*bt, *b, e);
                                if (BloqueSoportadoPoligono(nuevoB)) {
                                    opciones.push_back(nuevoB);
                                }
                            }
                        }
                    }
                    else if ((*bt).cantidad > 1) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                int desde = opciones.size();

                                // x

                                int nx = min(e.dx / (*b).x, (*bt).cantidad);
                                while (nx > 0) {
                                    Bloque nuevo(*bt, *b, e, nx, 1, 1);
                                    if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                        if (BloqueSoportadoPoligono(nuevo)) {
                                            opciones.push_back(nuevo);

                                            // y z

                                            int ny = min(e.dy / (*b).y, (*bt).cantidad / nx);
                                            while (ny > 1) {
                                                nuevo = Bloque(*bt, *b, e, nx, ny, 1);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        int nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                                        while (nz > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --nz;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --ny;
                                            }

                                            // z y

                                            int nz = min(e.dz / (*b).z, (*bt).cantidad / nx);
                                            while (nz > 1) {
                                                nuevo = Bloque(*bt, *b, e, nx, 1, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                                        while (ny > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --ny;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --nz;
                                            }

                                            break;
                                        }
                                    }
                                    --nx;
                                }

                                // y

                                int ny = min(e.dy / (*b).y, (*bt).cantidad);
                                while (ny > 0) {
                                    Bloque nuevo(*bt, *b, e, 1, ny, 1);
                                    if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                        if (BloqueSoportadoPoligono(nuevo)) {
                                            opciones.push_back(nuevo);

                                            // x z

                                            nx = min(e.dx / (*b).x, (*bt).cantidad / ny);
                                            while (nx > 1) {
                                                nuevo = Bloque(*bt, *b, e, nx, ny, 1);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        int nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                                        while (nz > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --nz;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --nx;
                                            }

                                            // z x

                                            int nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                            while (nz > 1) {
                                                nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        int nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                                        while (nx > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --nx;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --nz;
                                            }

                                            break;
                                        }
                                    }
                                    --ny;
                                }

                                // z

                                int nz = min(e.dz / (*b).z, (*bt).cantidad);
                                while (nz > 0) {
                                    Bloque nuevo(*bt, *b, e, 1, 1, nz);
                                    if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                        if (BloqueSoportadoPoligono(nuevo)) {
                                            opciones.push_back(nuevo);

                                            // x y

                                            nx = min(e.dx / (*b).x, (*bt).cantidad / nz);
                                            while (nx > 1) {
                                                nuevo = Bloque(*bt, *b, e, nx, 1, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        int ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                                        while (ny > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --ny;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --nx;
                                            }

                                            // y x

                                            int ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                            while (ny > 1) {
                                                nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                                        while (nx > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --nx;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --ny;
                                            }

                                            break;
                                        }
                                    }
                                    --nz;
                                }
                            }
                        }
                    }
                }
            }
        }
        else {
            if (dist01(generator) < 0.5) { // Bloques de 1 sola caja en x
                for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                    if ((*bt).cantidad == 1) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                opciones.push_back(Bloque(*bt, *b, e));
                            }
                        }
                    }
                    else if ((*bt).cantidad > 1) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                int desde = opciones.size();

                                // y z 

                                int ny = min(e.dy / (*b).y, (*bt).cantidad);
                                MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, 1));
                                int nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                if (nz > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, nz));
                                }

                                // z y

                                nz = min(e.dz / (*b).z, (*bt).cantidad);
                                MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, 1, nz));
                                ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                if (ny > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, nz));
                                }
                            }
                        }
                    }
                }
            }
            else { // Bloques de varias cajas en x
                for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                    if ((*bt).cantidad == 1) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                opciones.push_back(Bloque(*bt, *b, e));
                            }
                        }
                    }
                    else if ((*bt).cantidad > 1) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                int desde = opciones.size();

                                // x y z

                                int nx = min(e.dx / (*b).x, (*bt).cantidad);
                                MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, 1, 1));
                                int ny = min(e.dy / (*b).y, (*bt).cantidad / nx);
                                if (ny > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, 1));
                                    int nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                    if (nz > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }

                                // x z y

                                int nz = min(e.dz / (*b).z, (*bt).cantidad / nx);
                                if (nz > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, 1, nz));
                                    ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                    if (ny > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }

                                // y x z

                                ny = min(e.dy / (*b).y, (*bt).cantidad);
                                MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, 1));
                                nx = min(e.dx / (*b).x, (*bt).cantidad / ny);
                                if (nx > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, 1));
                                    nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                    if (nz > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }

                                // y z x

                                nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                if (nz > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, nz));
                                    nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                    if (nx > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }

                                // z x y

                                nz = min(e.dz / (*b).z, (*bt).cantidad);
                                MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, 1, nz));
                                nx = min(e.dx / (*b).x, (*bt).cantidad / nz);
                                if (nx > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, 1, nz));
                                    ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                    if (ny > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }

                                // z y x

                                ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                if (ny > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, nz));
                                    nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                    if (nx > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    void GenerarOpciones_Presion(vector<Bloque>& opciones, MaximalSpace const& e) {
        opciones.reserve(maxNBloques);
        if (dist01(generator) < 0.5) { // Bloques de 1 sola caja en x
            for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                if ((*bt).cantidad == 1) {
                    for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                        if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                            Bloque nuevoB(*bt, *b, e);
                            if (BloqueSoportadoPresion(nuevoB)) {
                                opciones.push_back(nuevoB);
                            }
                        }
                    }
                }
                else if ((*bt).cantidad > 1) {
                    for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                        if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                            int desde = opciones.size();

                            // y z

                            int ny = min(e.dy / (*b).y, (*bt).cantidad);
                            while (ny > 0) {
                                Bloque nuevo(*bt, *b, e, 1, ny, 1);
                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                    if (BloqueSoportadoPresion(nuevo)) {
                                        opciones.push_back(nuevo);
                                        int nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                        while (nz > 1) {
                                            nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    break;
                                                }
                                            }
                                            --nz;
                                        }
                                        break;
                                    }
                                }
                                --ny;
                            }

                            // z y

                            int nz = min(e.dz / (*b).z, (*bt).cantidad);
                            while (nz > 0) {
                                Bloque nuevo(*bt, *b, e, 1, 1, nz);
                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                    if (BloqueSoportadoPresion(nuevo)) {
                                        opciones.push_back(nuevo);
                                        int ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                        while (ny > 1) {
                                            nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    break;
                                                }
                                            }
                                            --ny;
                                        }
                                        break;
                                    }
                                }
                                --nz;
                            }
                        }
                    }
                }
            }
        }
        else { // Bloques de varias cajas en x
            for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                if ((*bt).cantidad == 1) {
                    for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                        if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                            Bloque nuevoB(*bt, *b, e);
                            if (BloqueSoportadoPresion(nuevoB)) {
                                opciones.push_back(nuevoB);
                            }
                        }
                    }
                }
                else if ((*bt).cantidad > 1) {
                    for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                        if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                            int desde = opciones.size();

                            // x

                            int nx = min(e.dx / (*b).x, (*bt).cantidad);
                            while (nx > 0) {
                                Bloque nuevo(*bt, *b, e, nx, 1, 1);
                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                    if (BloqueSoportadoPresion(nuevo)) {
                                        opciones.push_back(nuevo);

                                        // y z

                                        int ny = min(e.dy / (*b).y, (*bt).cantidad / nx);
                                        while (ny > 1) {
                                            nuevo = Bloque(*bt, *b, e, nx, ny, 1);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    int nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                                    while (nz > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --nz;
                                                    }
                                                    break;
                                                }
                                            }
                                            --ny;
                                        }

                                        // z y

                                        int nz = min(e.dz / (*b).z, (*bt).cantidad / nx);
                                        while (nz > 1) {
                                            nuevo = Bloque(*bt, *b, e, nx, 1, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                                    while (ny > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --ny;
                                                    }
                                                    break;
                                                }
                                            }
                                            --nz;
                                        }

                                        break;
                                    }
                                }
                                --nx;
                            }

                            // y

                            int ny = min(e.dy / (*b).y, (*bt).cantidad);
                            while (ny > 0) {
                                Bloque nuevo(*bt, *b, e, 1, ny, 1);
                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                    if (BloqueSoportadoPresion(nuevo)) {
                                        opciones.push_back(nuevo);

                                        // x z

                                        nx = min(e.dx / (*b).x, (*bt).cantidad / ny);
                                        while (nx > 1) {
                                            nuevo = Bloque(*bt, *b, e, nx, ny, 1);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    int nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                                    while (nz > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --nz;
                                                    }
                                                    break;
                                                }
                                            }
                                            --nx;
                                        }

                                        // z x

                                        int nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                        while (nz > 1) {
                                            nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    int nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                                    while (nx > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --nx;
                                                    }
                                                    break;
                                                }
                                            }
                                            --nz;
                                        }

                                        break;
                                    }
                                }
                                --ny;
                            }

                            // z

                            int nz = min(e.dz / (*b).z, (*bt).cantidad);
                            while (nz > 0) {
                                Bloque nuevo(*bt, *b, e, 1, 1, nz);
                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                    if (BloqueSoportadoPresion(nuevo)) {
                                        opciones.push_back(nuevo);

                                        // x y

                                        nx = min(e.dx / (*b).x, (*bt).cantidad / nz);
                                        while (nx > 1) {
                                            nuevo = Bloque(*bt, *b, e, nx, 1, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    int ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                                    while (ny > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --ny;
                                                    }
                                                    break;
                                                }
                                            }
                                            --nx;
                                        }

                                        // y x

                                        int ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                        while (ny > 1) {
                                            nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                                    while (nx > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --nx;
                                                    }
                                                    break;
                                                }
                                            }
                                            --ny;
                                        }

                                        break;
                                    }
                                }
                                --nz;
                            }
                        }
                    }
                }
            }
        }
    }
    void GenerarOpciones_Multidrop(vector<Bloque>& opciones, MaximalSpace const& e) {
        opciones.reserve(maxNBloques);
        if (r_estabilidad == 1) {
            if (dist01(generator) < 0.5) { // Bloques de 1 sola caja en x
                for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                    if ((*bt).cantidad == 1 && (*bt).cliente == clienteActual) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                Bloque nuevoB(*bt, *b, e);
                                if (BloqueSoportadoPoligono(nuevoB)) {
                                    opciones.push_back(nuevoB);
                                }
                            }
                        }
                    }
                    else if ((*bt).cantidad > 1 && (*bt).cliente == clienteActual) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                int desde = opciones.size();

                                // y z

                                int ny = min(e.dy / (*b).y, (*bt).cantidad);
                                while (ny > 0) {
                                    Bloque nuevo(*bt, *b, e, 1, ny, 1);
                                    if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                        if (BloqueSoportadoPoligono(nuevo)) {
                                            opciones.push_back(nuevo);
                                            int nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                            while (nz > 1) {
                                                nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        break;
                                                    }
                                                }
                                                --nz;
                                            }
                                            break;
                                        }
                                    }
                                    --ny;
                                }

                                // z y

                                int nz = min(e.dz / (*b).z, (*bt).cantidad);
                                while (nz > 0) {
                                    Bloque nuevo(*bt, *b, e, 1, 1, nz);
                                    if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                        if (BloqueSoportadoPoligono(nuevo)) {
                                            opciones.push_back(nuevo);
                                            int ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                            while (ny > 1) {
                                                nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        break;
                                                    }
                                                }
                                                --ny;
                                            }
                                            break;
                                        }
                                    }
                                    --nz;
                                }
                            }
                        }
                    }
                    else if ((*bt).cliente > clienteActual) {
                        break;
                    }
                }
            }
            else { // Bloques de varias cajas en x
                for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                    if ((*bt).cantidad == 1 && (*bt).cliente == clienteActual) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                Bloque nuevoB(*bt, *b, e);
                                if (BloqueSoportadoPoligono(nuevoB)) {
                                    opciones.push_back(nuevoB);
                                }
                            }
                        }
                    }
                    else if ((*bt).cantidad > 1 && (*bt).cliente == clienteActual) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                int desde = opciones.size();

                                // x

                                int nx = min(e.dx / (*b).x, (*bt).cantidad);
                                while (nx > 0) {
                                    Bloque nuevo(*bt, *b, e, nx, 1, 1);
                                    if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                        if (BloqueSoportadoPoligono(nuevo)) {
                                            opciones.push_back(nuevo);

                                            // y z

                                            int ny = min(e.dy / (*b).y, (*bt).cantidad / nx);
                                            while (ny > 1) {
                                                nuevo = Bloque(*bt, *b, e, nx, ny, 1);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        int nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                                        while (nz > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --nz;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --ny;
                                            }

                                            // z y

                                            int nz = min(e.dz / (*b).z, (*bt).cantidad / nx);
                                            while (nz > 1) {
                                                nuevo = Bloque(*bt, *b, e, nx, 1, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                                        while (ny > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --ny;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --nz;
                                            }

                                            break;
                                        }
                                    }
                                    --nx;
                                }

                                // y

                                int ny = min(e.dy / (*b).y, (*bt).cantidad);
                                while (ny > 0) {
                                    Bloque nuevo(*bt, *b, e, 1, ny, 1);
                                    if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                        if (BloqueSoportadoPoligono(nuevo)) {
                                            opciones.push_back(nuevo);

                                            // x z

                                            nx = min(e.dx / (*b).x, (*bt).cantidad / ny);
                                            while (nx > 1) {
                                                nuevo = Bloque(*bt, *b, e, nx, ny, 1);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        int nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                                        while (nz > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --nz;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --nx;
                                            }

                                            // z x

                                            int nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                            while (nz > 1) {
                                                nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        int nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                                        while (nx > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --nx;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --nz;
                                            }

                                            break;
                                        }
                                    }
                                    --ny;
                                }

                                // z

                                int nz = min(e.dz / (*b).z, (*bt).cantidad);
                                while (nz > 0) {
                                    Bloque nuevo(*bt, *b, e, 1, 1, nz);
                                    if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                        if (BloqueSoportadoPoligono(nuevo)) {
                                            opciones.push_back(nuevo);

                                            // x y

                                            nx = min(e.dx / (*b).x, (*bt).cantidad / nz);
                                            while (nx > 1) {
                                                nuevo = Bloque(*bt, *b, e, nx, 1, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        int ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                                        while (ny > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --ny;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --nx;
                                            }

                                            // y x

                                            int ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                            while (ny > 1) {
                                                nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                    if (BloqueSoportadoPoligono(nuevo)) {
                                                        opciones.push_back(nuevo);
                                                        nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                                        while (nx > 1) {
                                                            nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                                if (BloqueSoportadoPoligono(nuevo)) {
                                                                    opciones.push_back(nuevo);
                                                                    break;
                                                                }
                                                            }
                                                            --nx;
                                                        }
                                                        break;
                                                    }
                                                }
                                                --ny;
                                            }

                                            break;
                                        }
                                    }
                                    --nz;
                                }
                            }
                        }
                    }
                    else if ((*bt).cliente > clienteActual) {
                        break;
                    }
                }
            }
        }
        else {
            if (dist01(generator) < 0.5) { // Bloques de 1 sola caja en x
                for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                    if ((*bt).cantidad == 1 && (*bt).cliente == clienteActual) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                opciones.push_back(Bloque(*bt, *b, e));
                            }
                        }
                    }
                    else if ((*bt).cantidad > 1 && (*bt).cliente == clienteActual) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                int desde = opciones.size();

                                // y z 

                                int ny = min(e.dy / (*b).y, (*bt).cantidad);
                                MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, 1));
                                int nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                if (nz > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, nz));
                                }

                                // z y

                                nz = min(e.dz / (*b).z, (*bt).cantidad);
                                MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, 1, nz));
                                ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                if (ny > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, nz));
                                }
                            }
                        }
                    }
                    else if ((*bt).cliente > clienteActual) {
                        break;
                    }
                }
            }
            else { // Bloques de varias cajas en x
                for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                    if ((*bt).cantidad == 1 && (*bt).cliente == clienteActual) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                opciones.push_back(Bloque(*bt, *b, e));
                            }
                        }
                    }
                    else if ((*bt).cantidad > 1 && (*bt).cliente == clienteActual) {
                        for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                            if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                                int desde = opciones.size();

                                // x y z

                                int nx = min(e.dx / (*b).x, (*bt).cantidad);
                                MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, 1, 1));
                                int ny = min(e.dy / (*b).y, (*bt).cantidad / nx);
                                if (ny > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, 1));
                                    int nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                    if (nz > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }

                                // x z y

                                int nz = min(e.dz / (*b).z, (*bt).cantidad / nx);
                                if (nz > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, 1, nz));
                                    ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                    if (ny > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }

                                // y x z

                                ny = min(e.dy / (*b).y, (*bt).cantidad);
                                MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, 1));
                                nx = min(e.dx / (*b).x, (*bt).cantidad / ny);
                                if (nx > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, 1));
                                    nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                    if (nz > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }

                                // y z x

                                nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                if (nz > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, nz));
                                    nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                    if (nx > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }

                                // z x y

                                nz = min(e.dz / (*b).z, (*bt).cantidad);
                                MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, 1, nz));
                                nx = min(e.dx / (*b).x, (*bt).cantidad / nz);
                                if (nx > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, 1, nz));
                                    ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                    if (ny > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }

                                // z y x

                                ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                if (ny > 1) {
                                    MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, 1, ny, nz));
                                    nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                    if (nx > 1) {
                                        MeterBloqueNuevo(opciones, desde, Bloque(*bt, *b, e, nx, ny, nz));
                                    }
                                }
                            }
                        }
                    }
                    else if ((*bt).cliente > clienteActual) {
                        break;
                    }
                }
            }
        }
    }
    void GenerarOpciones_Multidrop_Presion(vector<Bloque>& opciones, MaximalSpace const& e) {
        opciones.reserve(maxNBloques);
        if (dist01(generator) < 0.5) { // Bloques de 1 sola caja en x
            for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                if ((*bt).cantidad == 1 && (*bt).cliente == clienteActual) {
                    for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                        if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                            Bloque nuevoB(*bt, *b, e);
                            if (BloqueSoportadoPresion(nuevoB)) {
                                opciones.push_back(nuevoB);
                            }
                        }
                    }
                }
                else if ((*bt).cantidad > 1 && (*bt).cliente == clienteActual) {
                    for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                        if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                            int desde = opciones.size();

                            // y z

                            int ny = min(e.dy / (*b).y, (*bt).cantidad);
                            while (ny > 0) {
                                Bloque nuevo(*bt, *b, e, 1, ny, 1);
                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                    if (BloqueSoportadoPresion(nuevo)) {
                                        opciones.push_back(nuevo);
                                        int nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                        while (nz > 1) {
                                            nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    break;
                                                }
                                            }
                                            --nz;
                                        }
                                        break;
                                    }
                                }
                                --ny;
                            }

                            // z y

                            int nz = min(e.dz / (*b).z, (*bt).cantidad);
                            while (nz > 0) {
                                Bloque nuevo(*bt, *b, e, 1, 1, nz);
                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                    if (BloqueSoportadoPresion(nuevo)) {
                                        opciones.push_back(nuevo);
                                        int ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                        while (ny > 1) {
                                            nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    break;
                                                }
                                            }
                                            --ny;
                                        }
                                        break;
                                    }
                                }
                                --nz;
                            }
                        }
                    }
                }
                else if ((*bt).cliente > clienteActual) {
                    break;
                }
            }
        }
        else { // Bloques de varias cajas en x
            for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
                if ((*bt).cantidad == 1 && (*bt).cliente == clienteActual) {
                    for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                        if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                            Bloque nuevoB(*bt, *b, e);
                            if (BloqueSoportadoPresion(nuevoB)) {
                                opciones.push_back(nuevoB);
                            }
                        }
                    }
                }
                else if ((*bt).cantidad > 1 && (*bt).cliente == clienteActual) {
                    for (vector<Box>::iterator b = (*bt).boxes.begin(); b < (*bt).boxes.end(); ++b) {
                        if ((*b).x <= e.dx && (*b).y <= e.dy && (*b).z <= e.dz) {
                            int desde = opciones.size();

                            // x

                            int nx = min(e.dx / (*b).x, (*bt).cantidad);
                            while (nx > 0) {
                                Bloque nuevo(*bt, *b, e, nx, 1, 1);
                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                    if (BloqueSoportadoPresion(nuevo)) {
                                        opciones.push_back(nuevo);

                                        // y z

                                        int ny = min(e.dy / (*b).y, (*bt).cantidad / nx);
                                        while (ny > 1) {
                                            nuevo = Bloque(*bt, *b, e, nx, ny, 1);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    int nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                                    while (nz > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --nz;
                                                    }
                                                    break;
                                                }
                                            }
                                            --ny;
                                        }

                                        // z y

                                        int nz = min(e.dz / (*b).z, (*bt).cantidad / nx);
                                        while (nz > 1) {
                                            nuevo = Bloque(*bt, *b, e, nx, 1, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                                    while (ny > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --ny;
                                                    }
                                                    break;
                                                }
                                            }
                                            --nz;
                                        }
                                        break;
                                    }
                                }
                                --nx;
                            }

                            // y

                            int ny = min(e.dy / (*b).y, (*bt).cantidad);
                            while (ny > 0) {
                                Bloque nuevo(*bt, *b, e, 1, ny, 1);
                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                    if (BloqueSoportadoPresion(nuevo)) {
                                        opciones.push_back(nuevo);

                                        // x z

                                        nx = min(e.dx / (*b).x, (*bt).cantidad / ny);
                                        while (nx > 1) {
                                            nuevo = Bloque(*bt, *b, e, nx, ny, 1);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    int nz = min(e.dz / (*b).z, (*bt).cantidad / (nx * ny));
                                                    while (nz > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --nz;
                                                    }
                                                    break;
                                                }
                                            }
                                            --nx;
                                        }

                                        // z x

                                        int nz = min(e.dz / (*b).z, (*bt).cantidad / ny);
                                        while (nz > 1) {
                                            nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    int nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                                    while (nx > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --nx;
                                                    }
                                                    break;
                                                }
                                            }
                                            --nz;
                                        }

                                        break;
                                    }
                                }
                                --ny;
                            }

                            // z

                            int nz = min(e.dz / (*b).z, (*bt).cantidad);
                            while (nz > 0) {
                                Bloque nuevo(*bt, *b, e, 1, 1, nz);
                                if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                    if (BloqueSoportadoPresion(nuevo)) {
                                        opciones.push_back(nuevo);

                                        // x y

                                        nx = min(e.dx / (*b).x, (*bt).cantidad / nz);
                                        while (nx > 1) {
                                            nuevo = Bloque(*bt, *b, e, nx, 1, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    int ny = min(e.dy / (*b).y, (*bt).cantidad / (nx * nz));
                                                    while (ny > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --ny;
                                                    }
                                                    break;
                                                }
                                            }
                                            --nx;
                                        }

                                        // y x

                                        int ny = min(e.dy / (*b).y, (*bt).cantidad / nz);
                                        while (ny > 1) {
                                            nuevo = Bloque(*bt, *b, e, 1, ny, nz);
                                            if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                if (BloqueSoportadoPresion(nuevo)) {
                                                    opciones.push_back(nuevo);
                                                    nx = min(e.dx / (*b).x, (*bt).cantidad / (ny * nz));
                                                    while (nx > 1) {
                                                        nuevo = Bloque(*bt, *b, e, nx, ny, nz);
                                                        if (EsBloqueNuevo(opciones, desde, nuevo)) {
                                                            if (BloqueSoportadoPresion(nuevo)) {
                                                                opciones.push_back(nuevo);
                                                                break;
                                                            }
                                                        }
                                                        --nx;
                                                    }
                                                    break;
                                                }
                                            }
                                            --ny;
                                        }

                                        break;
                                    }
                                }
                                --nz;
                            }
                        }
                    }
                }
                else if ((*bt).cliente > clienteActual) {
                    break;
                }
            }
        }
    }
    void SeleccionarBloque_BestFit_Deterministico(vector<Bloque>& opciones, vector<Bloque>::iterator& resp) {
        for (vector<Bloque>::iterator b_it = opciones.begin() + 1; b_it < opciones.end(); ++b_it) {
            if ((*resp).bf1 < (*b_it).bf1) resp = b_it;
            else if ((*resp).bf1 == (*b_it).bf1) {
                if ((*resp).bf2 < (*b_it).bf2) resp = b_it;
                else if ((*resp).bf2 == (*b_it).bf2) {
                    if ((*resp).bf3 < (*b_it).bf3) resp = b_it;
                    else if ((*resp).bf3 == (*b_it).bf3) {
                        if ((*resp).vol < (*b_it).vol) resp = b_it;
                    }
                }
            }
        }
    }
    void SeleccionarBloque_BestFit(vector<Bloque>& opciones, vector<Bloque>::iterator& resp) {
        sort(opciones.begin(), opciones.end(), Orden_BestFit_MaxVol);
        uniform_int_distribution<int> distri(0, miBusquedaBinariaBF(opciones, opciones.front().bf1 + (double)(opciones.back().bf1 - opciones.front().bf1) * alpha + error));
        resp = opciones.begin() + distri(generator);
    }
    void SeleccionarBloque_MaxVol(vector<Bloque>& opciones, vector<Bloque>::iterator& resp) {
        sort(opciones.begin(), opciones.end(), Orden_MaxVol_BestFit);
        uniform_int_distribution<int> distri(0, miBusquedaBinariaVol(opciones, opciones.front().vol + (double)(opciones.back().vol - opciones.front().vol) * alpha + error));
        resp = opciones.begin() + distri(generator);
    }

    // Empacar

    void ActualizarArbolPresion(int const& nUltimasPiezas) {

        // Determinar arbol de actualización

        vector<vector<PackedBox>::iterator> arbol;
        int menorZ1 = ContenedorDimz;
        for (vector<PackedBox>::iterator p_it = empacados.begin() + (int)empacados.size() - nUltimasPiezas; p_it < empacados.end(); ++p_it) {
            arbol.push_back(p_it);
            if ((*p_it).z1 < menorZ1) {
                menorZ1 = (*p_it).z1;
            }
        }
        int hasta = arbol.size();
        for (int i = 0; i < hasta; ++i) {
            vector<PackedBox>::iterator p_it = arbol[i];
            if ((*p_it).z1 == menorZ1) {
                int siguiente = arbol.size();
                while ((*p_it).z1 > 0) {
                    for (vector<Soporte>::iterator s_it = (*p_it).soportes.begin(); s_it < (*p_it).soportes.end(); ++s_it) {
                        bool nuevaPieza = true;
                        for (vector<vector<PackedBox>::iterator>::iterator p_itit = arbol.begin() + hasta; p_itit < arbol.end(); ++p_itit) {
                            if ((*s_it).indPB == (**p_itit).ind) {
                                nuevaPieza = false;
                                break;
                            }
                        }
                        if (nuevaPieza) {
                            arbol.push_back(empacados.begin() + (*s_it).indPB);
                        }
                    }
                    if (siguiente < arbol.size()) {
                        for (; siguiente < arbol.size(); ++siguiente) {
                            p_it = arbol[siguiente];
                            if ((*p_it).z1 > 0) break;
                        }
                        ++siguiente;
                    }
                    else break;
                }
            }
        }
        sort(arbol.begin(), arbol.end(), [](vector<PackedBox>::iterator const& p1, vector<PackedBox>::iterator const& p2)->bool {return (*p1).z1 > (*p2).z1; });

        // Actualizar soporte

        for (vector<vector<PackedBox>::iterator>::iterator p_itit = arbol.begin(); p_itit < arbol.end(); ++p_itit) {
            if ((**p_itit).z1 > 0) {
                for (vector<Soporte>::iterator s_it = (**p_itit).soportes.begin(); s_it < (**p_itit).soportes.end(); ++s_it) {

                    // Actualizar soporte

                    vector<vector<PackedBox>::iterator>::iterator p_itit2 = p_itit + 1;
                    for (; p_itit2 < arbol.end(); ++p_itit2) {
                        if ((**p_itit2).ind == (*s_it).indPB) {
                            break;
                        }
                    }
                    (**p_itit2).ActualizarCentroDeMasa(*s_it);
                }
            }
        }

        // Actualizar distribución de pesos

        for (vector<vector<PackedBox>::iterator>::iterator p_itit = arbol.begin(); p_itit < arbol.end(); ++p_itit) {
            if ((**p_itit).z1 > 0) {
                (**p_itit).ActualizarSoportes();
            }
        }
    }
    void Empacar(vector<Bloque>::iterator& b, vector<MaximalSpace>::iterator& e) {
        for (vector<BloqueAux>::iterator b_it = (*b).cajas.begin(); b_it < (*b).cajas.end(); ++b_it) {
            empacados.push_back(PackedBox(empacados.size(), *b, *b_it));
        }
        utilizacion += (*b).bt.volumen * (*b).btq;
        boxest[(*b).bt.id].cantidad -= (*b).btq;
        if (boxest[(*b).bt.id].cantidad == 0) { // Se revisa si hay cajas por empacar
            hayCajasPorEmpacar = false;
            for (vector<BoxType>::iterator b = boxest.begin(); b < boxest.end(); ++b) {
                if ((*b).cantidad > 0) {
                    hayCajasPorEmpacar = true;
                    break;
                }
            }
            if (!hayCajasPorEmpacar) return;
            if ((*b).bt.minDx == minDx || (*b).bt.minDy == minDy || (*b).bt.minDz == minDz) ActualizarMinD();
        }
        JuntarEspacios(ActualizarEspacios(*b));
        EliminarEspaciosMinD();
        if (r_estabilidad == 1) {
            ActualizarIndicesCajasDeEspacios();
        }
    }
    void Empacar_Multidrop(vector<Bloque>::iterator& b, vector<MaximalSpace>::iterator& e) {
        for (vector<BloqueAux>::iterator b_it = (*b).cajas.begin(); b_it < (*b).cajas.end(); ++b_it) {
            empacados.push_back(PackedBox(empacados.size(), *b, *b_it));
        }
        utilizacion += (*b).bt.volumen * (*b).btq;
        boxest[(*b).bt.id].cantidad -= (*b).btq;
        JuntarEspacios(ActualizarEspacios(*b));
        if (boxest[(*b).bt.id].cantidad == 0) {

            // Se revisa si hay cajas por empacar

            hayCajasPorEmpacar = false;
            for (vector<BoxType>::iterator b = boxest.begin(); b < boxest.end(); ++b) {
                if ((*b).cantidad > 0) {
                    hayCajasPorEmpacar = true;
                    break;
                }
            }
            if (!hayCajasPorEmpacar) return;

            // Se analiza si el cliente actual cambió

            bool hayMasCajasDelCliente = false;
            for (vector<BoxType>::iterator bt_it = boxest.begin() + (*b).bt.id; bt_it < boxest.end(); ++bt_it) {
                if ((*bt_it).cliente != clienteActual) break;
                else if ((*bt_it).cantidad > 0) {
                    hayMasCajasDelCliente = true;
                    break;
                }
            }
            if (!hayMasCajasDelCliente) {
                for (vector<BoxType>::reverse_iterator bt_rit = boxest.rbegin() + (int)boxest.size() - (*b).bt.id; bt_rit != boxest.rend(); ++bt_rit) {
                    if ((*bt_rit).cliente != clienteActual) break;
                    else if ((*bt_rit).cantidad > 0) {
                        hayMasCajasDelCliente = true;
                        break;
                    }
                }
            }
            if (!hayMasCajasDelCliente) {
                ++clienteActual;
                if (r_multidrop == 1) ActualizarEspaciosVisibilidad();
                else if (r_multidrop == 2) {
                    int x0 = empacados.back().x2;
                    for (vector<PackedBox>::reverse_iterator p = empacados.rbegin() + 1; p != empacados.rend(); ++p) {
                        if ((*p).x2 > x0) x0 = (*p).x2;
                    }
                    ActualizarEspaciosCapas(x0);
                }
            }
            if ((*b).bt.minDx == minDx || (*b).bt.minDy == minDy || (*b).bt.minDz == minDz) ActualizarMinD();
        }
        EliminarEspaciosMinD();
        if (r_estabilidad == 1) {
            ActualizarIndicesCajasDeEspacios();
        }
    }
    void GenerarPiezasEmpacadas(Bloque& b) {

        // Generar las piezas empacadas

        int desde = empacados.size();
        for (vector<BloqueAux>::iterator b_it = b.cajas.begin(); b_it < b.cajas.end(); ++b_it) {
            if ((*b_it).z1 == b.z1) {
                if (b.z1 == 0) {
                    empacados.push_back(PackedBox(empacados.size(), b, *b_it));
                }
                else {
                    double centroMasax = (double)((*b_it).x1 + (double)b.b.x / 2.0);
                    double centroMasay = (double)((*b_it).y1 + (double)b.b.y / 2.0);
                    vector<Point> intersecciones;
                    intersecciones.reserve(4 * b.indPB_EM.size());
                    vector<Soporte> misSoportes; misSoportes.reserve(b.indPB_EM.size());
                    for (vector<int>::iterator idPB = b.indPB_EM.begin(); idPB < b.indPB_EM.end(); ++idPB) {
                        PackedBox& p = empacados[*idPB];
                        if (!((*b_it).x1 >= p.x2 || p.x1 >= (*b_it).x1 + b.b.x)) {
                            if (!((*b_it).y1 >= p.y2 || p.y1 >= (*b_it).y1 + b.b.y)) {

                                // Se encuentran los puntos de intersección

                                vector<int> misx({ (*b_it).x1, (*b_it).x1 + b.b.x, p.x1, p.x2 });
                                sort(misx.begin(), misx.end());
                                vector<int> misy({ (*b_it).y1, (*b_it).y1 + b.b.y, p.y1, p.y2 });
                                sort(misy.begin(), misy.end());
                                intersecciones.push_back(Point(misx[1], misy[1]));
                                intersecciones.push_back(Point(misx[2], misy[1]));
                                intersecciones.push_back(Point(misx[1], misy[2]));
                                intersecciones.push_back(Point(misx[2], misy[2]));

                                // Se crea el soporte

                                misSoportes.push_back(Soporte(*idPB, (misx[2] - misx[1]) * (misy[2] - misy[1]), (double)(misx[2] + misx[1]) / 2.0f, (double)(misy[2] + misy[1]) / 2.0f, p.presion));
                            }
                        }
                    }

                    // Determinar polígono

                    if (intersecciones.size() > 0) {
                        Poligono poligono(intersecciones);
                        empacados.push_back(PackedBox(empacados.size(), b, *b_it, centroMasax, centroMasay, poligono, misSoportes));
                        if (empacados.back().z2 == b.z2 && empacados.back().z1 > 0) {
                            empacados.back().ActualizarSoportes();
                        }
                    }
                    else cout << "Error" << endl;
                }
            }
            else {
                double centroMasax = (double)((*b_it).x1 + (double)b.b.x / 2.0);
                double centroMasay = (double)((*b_it).y1 + (double)b.b.y / 2.0);
                vector<Point> intersecciones;
                intersecciones.reserve(4 * b.indPB_EM.size());
                vector<Soporte> misSoportes; misSoportes.reserve(b.indPB_EM.size());
                for (vector<PackedBox>::iterator p_it = empacados.begin() + desde; p_it < empacados.end(); ++p_it) {
                    if (!((*b_it).x1 >= (*p_it).x2 || (*p_it).x1 >= (*b_it).x1 + b.b.x)) {
                        if (!((*b_it).y1 >= (*p_it).y2 || (*p_it).y1 >= (*b_it).y1 + b.b.y)) {

                            // Se encuentran los puntos de intersección

                            vector<int> misx({ (*b_it).x1, (*b_it).x1 + b.b.x, (*p_it).x1, (*p_it).x2 });
                            sort(misx.begin(), misx.end());
                            vector<int> misy({ (*b_it).y1, (*b_it).y1 + b.b.y, (*p_it).y1, (*p_it).y2 });
                            sort(misy.begin(), misy.end());
                            intersecciones.push_back(Point(misx[1], misy[1]));
                            intersecciones.push_back(Point(misx[2], misy[1]));
                            intersecciones.push_back(Point(misx[1], misy[2]));
                            intersecciones.push_back(Point(misx[2], misy[2]));

                            // Se crea el soporte

                            misSoportes.push_back(Soporte((*p_it).ind, (misx[2] - misx[1]) * (misy[2] - misy[1]), (double)(misx[2] + misx[1]) / 2.0f, (double)(misy[2] + misy[1]) / 2.0f, (*p_it).presion));
                        }
                    }
                }

                // Determinar polígono

                if (intersecciones.size() > 0) {
                    Poligono poligono(intersecciones);
                    empacados.push_back(PackedBox(empacados.size(), b, *b_it, centroMasax, centroMasay, poligono, misSoportes));
                    if (empacados.back().z2 == b.z2 && empacados.back().z1 > 0) {
                        empacados.back().ActualizarSoportes();
                    }
                }
                else cout << "Error" << endl;
            }
        }
    }
    void Empacar_Presion(vector<Bloque>::iterator& b, vector<MaximalSpace>::iterator& e) {
        GenerarPiezasEmpacadas(*b);
        ActualizarArbolPresion((*b).btq);
        utilizacion += (*b).bt.volumen * (*b).btq;
        boxest[(*b).bt.id].cantidad -= (*b).btq;
        if (boxest[(*b).bt.id].cantidad == 0) { // Se revisa si hay cajas por empacar
            hayCajasPorEmpacar = false;
            for (vector<BoxType>::iterator b = boxest.begin(); b < boxest.end(); ++b) {
                if ((*b).cantidad > 0) {
                    hayCajasPorEmpacar = true;
                    break;
                }
            }
            if (!hayCajasPorEmpacar) return;
            if ((*b).bt.minDx == minDx || (*b).bt.minDy == minDy || (*b).bt.minDz == minDz) ActualizarMinD();
        }
        JuntarEspacios(ActualizarEspacios(*b));
        EliminarEspaciosMinD();
        ActualizarIndicesCajasDeEspacios();
    }
    void Empacar_Multidrop_Presion(vector<Bloque>::iterator& b, vector<MaximalSpace>::iterator& e) {
        GenerarPiezasEmpacadas(*b);
        ActualizarArbolPresion((*b).cajas.size());
        utilizacion += (*b).bt.volumen * (*b).btq;
        boxest[(*b).bt.id].cantidad -= (*b).btq;
        JuntarEspacios(ActualizarEspacios(*b));
        if (boxest[(*b).bt.id].cantidad == 0) {

            // Se revisa si hay cajas por empacar

            hayCajasPorEmpacar = false;
            for (vector<BoxType>::iterator b = boxest.begin(); b < boxest.end(); ++b) {
                if ((*b).cantidad > 0) {
                    hayCajasPorEmpacar = true;
                    break;
                }
            }
            if (!hayCajasPorEmpacar) return;

            // Se analiza si el cliente actual cambió

            bool hayMasCajasDelCliente = false;
            for (vector<BoxType>::iterator bt_it = boxest.begin() + (*b).bt.id; bt_it < boxest.end(); ++bt_it) {
                if ((*bt_it).cliente != clienteActual) break;
                else if ((*bt_it).cantidad > 0) {
                    hayMasCajasDelCliente = true;
                    break;
                }
            }
            if (!hayMasCajasDelCliente) {
                for (vector<BoxType>::reverse_iterator bt_rit = boxest.rbegin() + (int)boxest.size() - (*b).bt.id; bt_rit != boxest.rend(); ++bt_rit) {
                    if ((*bt_rit).cliente != clienteActual) break;
                    else if ((*bt_rit).cantidad > 0) {
                        hayMasCajasDelCliente = true;
                        break;
                    }
                }
            }
            if (!hayMasCajasDelCliente) {
                ++clienteActual;
                if (r_multidrop == 1) ActualizarEspaciosVisibilidad();
                else if (r_multidrop == 2) {
                    int x0 = empacados.back().x2;
                    for (vector<PackedBox>::reverse_iterator p = empacados.rbegin() + 1; p != empacados.rend(); ++p) {
                        if ((*p).x2 > x0) x0 = (*p).x2;
                    }
                    ActualizarEspaciosCapas(x0);
                }
            }
            if ((*b).bt.minDx == minDx || (*b).bt.minDy == minDy || (*b).bt.minDz == minDz) ActualizarMinD();
        }
        EliminarEspaciosMinD();
        ActualizarIndicesCajasDeEspacios();
        if (holi == 16) {
            hayCajasPorEmpacar = false;
        }
    }

    // Constructivo

public:

    void Constructivo() {
        //++holi;
        //cout << holi << endl;
        while ((int)espacios.size() > 0 && hayCajasPorEmpacar) {

            // Se selecciona el espacio maximal

            vector<MaximalSpace>::iterator mejorEspacio = espacios.begin();
            double myRand1 = dist01(generator);
            if (myRand1 < 0.33) SeleccionarEspacio_Fondo_Altura(mejorEspacio);
            else if (myRand1 < 0.67) SeleccionarEspacio_Fondo_Bajo(mejorEspacio);
            else SeleccionarEspacio_Altura_Fondo(mejorEspacio);

            // Se selecciona la pieza a empacar

            vector<Bloque> opciones;
            GenerarOpciones(opciones, *mejorEspacio);
            if (opciones.size() > 0) {
                vector<Bloque>::iterator mejorBloque = opciones.begin();
                if (dist01(generator) < 0.5) SeleccionarBloque_BestFit(opciones, mejorBloque);
                else SeleccionarBloque_MaxVol(opciones, mejorBloque);
                Empacar(mejorBloque, mejorEspacio);
            }
            else espacios.erase(mejorEspacio);
        }
    }
    void ConstructivoDeterministico() {
        while ((int)espacios.size() > 0 && hayCajasPorEmpacar) {

            // Se selecciona el espacio maximal

            vector<MaximalSpace>::iterator mejorEspacio = espacios.begin();
            double myRand1 = dist01(generator);
            if (myRand1 < 0.33) SeleccionarEspacio_Fondo_Altura(mejorEspacio);
            else if (myRand1 < 0.67) SeleccionarEspacio_Fondo_Bajo(mejorEspacio);
            else SeleccionarEspacio_Altura_Fondo(mejorEspacio);

            // Se selecciona la pieza a empacar

            vector<Bloque> opciones;
            GenerarOpciones(opciones, *mejorEspacio);
            if (opciones.size() > 0) {
                vector<Bloque>::iterator mejorBloque = opciones.begin();
                SeleccionarBloque_BestFit_Deterministico(opciones, mejorBloque);
                Empacar(mejorBloque, mejorEspacio);
            }
            else espacios.erase(mejorEspacio);
        }
    }
    void Constructivo_Multidrop() {
        //++holi;
        while ((int)espacios.size() > 0 && hayCajasPorEmpacar) { // Falta poner restricción si hay algún tipo con caja

            // Se selecciona el espacio maximal

            vector<MaximalSpace>::iterator mejorEspacio = espacios.begin();
            double myRand1 = dist01(generator);
            if (myRand1 < 0.33) SeleccionarEspacio_Fondo_Altura(mejorEspacio);
            else if (myRand1 < 0.67) SeleccionarEspacio_Fondo_Bajo(mejorEspacio);
            else SeleccionarEspacio_Altura_Fondo(mejorEspacio);

            // Se selecciona la pieza a empacar

            vector<Bloque> opciones;
            GenerarOpciones_Multidrop(opciones, *mejorEspacio);
            double myRand2 = (double)dist01(generator);
            if (opciones.size() > 0) {
                vector<Bloque>::iterator mejorBloque = opciones.begin();
                if (myRand2 < 0.5) SeleccionarBloque_BestFit(opciones, mejorBloque);
                else SeleccionarBloque_MaxVol(opciones, mejorBloque);
                Empacar_Multidrop(mejorBloque, mejorEspacio);
            }
            else {
                if (myRand1 < 0.33) sort(espacios.begin(), espacios.end(), Orden_Fondo_Altura);
                else if (myRand1 < 0.67) sort(espacios.begin(), espacios.end(), Orden_Fondo_Bajo);
                else sort(espacios.begin(), espacios.end(), Orden_Altura_Fondo);

                // Recorrer los demás espacios en orden

                if (myRand2 < 0.5f) {
                    bool terminar = true;
                    for (vector<MaximalSpace>::iterator e = espacios.begin() + 1; e < espacios.end(); ++e) {
                        GenerarOpciones_Multidrop(opciones, *e);
                        if (opciones.size() > 0) {
                            vector<Bloque>::iterator mejorBloque = opciones.begin();
                            SeleccionarBloque_BestFit(opciones, mejorBloque);
                            Empacar_Multidrop(mejorBloque, e);
                            terminar = false;
                            break;
                        }
                    }
                    if (terminar) break;
                }
                else {
                    bool terminar = true;
                    for (vector<MaximalSpace>::iterator e = espacios.begin() + 1; e < espacios.end(); ++e) {
                        GenerarOpciones_Multidrop(opciones, *e);
                        if (opciones.size() > 0) {
                            vector<Bloque>::iterator mejorBloque = opciones.begin();
                            SeleccionarBloque_MaxVol(opciones, mejorBloque);
                            Empacar_Multidrop(mejorBloque, e);
                            terminar = false;
                            break;
                        }
                    }
                    if (terminar) break;
                }
            }
        }
    }
    void ConstructivoDeterministico_Multidrop() {
        while ((int)espacios.size() > 0 && hayCajasPorEmpacar) {

            // Se selecciona el espacio maximal

            vector<MaximalSpace>::iterator mejorEspacio = espacios.begin();
            double myRand1 = dist01(generator);
            if (myRand1 < 0.33) SeleccionarEspacio_Fondo_Altura(mejorEspacio);
            else if (myRand1 < 0.67) SeleccionarEspacio_Fondo_Bajo(mejorEspacio);
            else SeleccionarEspacio_Altura_Fondo(mejorEspacio);

            // Se selecciona la pieza a empacar

            vector<Bloque> opciones;
            GenerarOpciones_Multidrop(opciones, *mejorEspacio);
            double myRand2 = (double)dist01(generator);
            if (opciones.size() > 0) {
                vector<Bloque>::iterator mejorBloque = opciones.begin();
                SeleccionarBloque_BestFit_Deterministico(opciones, mejorBloque);
                Empacar_Multidrop(mejorBloque, mejorEspacio);
            }
            else {
                if (myRand1 < 0.33) sort(espacios.begin(), espacios.end(), Orden_Fondo_Altura);
                else if (myRand1 < 0.67) sort(espacios.begin(), espacios.end(), Orden_Fondo_Bajo);
                else sort(espacios.begin(), espacios.end(), Orden_Altura_Fondo);

                // Recorrer los demás espacios en orden

                bool terminar = true;
                for (vector<MaximalSpace>::iterator e = espacios.begin() + 1; e < espacios.end(); ++e) {
                    GenerarOpciones_Multidrop(opciones, *e);
                    if (opciones.size() > 0) {
                        vector<Bloque>::iterator mejorBloque = opciones.begin();
                        SeleccionarBloque_BestFit_Deterministico(opciones, mejorBloque);
                        Empacar_Multidrop(mejorBloque, e);
                        terminar = false;
                        break;
                    }
                }
                if (terminar) break;
            }
        }
    }
    void Constructivo_Presion() {
        //++holi;
        while ((int)espacios.size() > 0 && hayCajasPorEmpacar) {

            // Se selecciona el espacio maximal

            vector<MaximalSpace>::iterator mejorEspacio = espacios.begin();
            double myRand1 = dist01(generator);
            if (myRand1 < 0.33) SeleccionarEspacio_Fondo_Altura(mejorEspacio);
            else if (myRand1 < 0.67) SeleccionarEspacio_Fondo_Bajo(mejorEspacio);
            else SeleccionarEspacio_Altura_Fondo(mejorEspacio);

            // Se selecciona la pieza a empacar

            vector<Bloque> opciones;
            GenerarOpciones_Presion(opciones, *mejorEspacio);
            if (opciones.size() > 0) {
                vector<Bloque>::iterator mejorBloque = opciones.begin();
                if (dist01(generator) < 0.5) SeleccionarBloque_BestFit(opciones, mejorBloque);
                else SeleccionarBloque_MaxVol(opciones, mejorBloque);
                Empacar_Presion(mejorBloque, mejorEspacio);
            }
            else espacios.erase(mejorEspacio);
        }
    }
    void ConstructivoDeterministico_Presion() {
        while ((int)espacios.size() > 0 && hayCajasPorEmpacar) {

            // Se selecciona el espacio maximal

            vector<MaximalSpace>::iterator mejorEspacio = espacios.begin();
            double myRand1 = dist01(generator);
            if (myRand1 < 0.33) SeleccionarEspacio_Fondo_Altura(mejorEspacio);
            else if (myRand1 < 0.67) SeleccionarEspacio_Fondo_Bajo(mejorEspacio);
            else SeleccionarEspacio_Altura_Fondo(mejorEspacio);

            // Se selecciona la pieza a empacar

            vector<Bloque> opciones;
            GenerarOpciones_Presion(opciones, *mejorEspacio);
            if (opciones.size() > 0) {
                vector<Bloque>::iterator mejorBloque = opciones.begin();
                SeleccionarBloque_BestFit_Deterministico(opciones, mejorBloque);
                Empacar_Presion(mejorBloque, mejorEspacio);
            }
            else espacios.erase(mejorEspacio);
        }
    }
    void Constructivo_Multidrop_Presion() {
        //++holi;
        while ((int)espacios.size() > 0 && hayCajasPorEmpacar) { // Falta poner restricción si hay algún tipo con caja

            // Se selecciona el espacio maximal

            vector<MaximalSpace>::iterator mejorEspacio = espacios.begin();
            double myRand1 = dist01(generator);
            if (myRand1 < 0.33) SeleccionarEspacio_Fondo_Altura(mejorEspacio);
            else if (myRand1 < 0.67) SeleccionarEspacio_Fondo_Bajo(mejorEspacio);
            else SeleccionarEspacio_Altura_Fondo(mejorEspacio);

            // Se selecciona la pieza a empacar

            vector<Bloque> opciones;
            GenerarOpciones_Multidrop_Presion(opciones, *mejorEspacio);
            double myRand2 = (double)dist01(generator);
            if (opciones.size() > 0) {
                vector<Bloque>::iterator mejorBloque = opciones.begin();
                if (myRand2 < 0.5) SeleccionarBloque_BestFit(opciones, mejorBloque);
                else SeleccionarBloque_MaxVol(opciones, mejorBloque);
                Empacar_Multidrop_Presion(mejorBloque, mejorEspacio);
            }
            else {
                if (myRand1 < 0.33) sort(espacios.begin(), espacios.end(), Orden_Fondo_Altura);
                else if (myRand1 < 0.67) sort(espacios.begin(), espacios.end(), Orden_Fondo_Bajo);
                else sort(espacios.begin(), espacios.end(), Orden_Altura_Fondo);

                // Recorrer los demás espacios en orden

                if (myRand2 < 0.5f) {
                    bool terminar = true;
                    for (vector<MaximalSpace>::iterator e = espacios.begin() + 1; e < espacios.end(); ++e) {
                        GenerarOpciones_Multidrop_Presion(opciones, *e);
                        if (opciones.size() > 0) {
                            vector<Bloque>::iterator mejorBloque = opciones.begin();
                            SeleccionarBloque_BestFit(opciones, mejorBloque);
                            Empacar_Multidrop_Presion(mejorBloque, e);
                            terminar = false;
                            break;
                        }
                    }
                    if (terminar) break;
                }
                else {
                    bool terminar = true;
                    for (vector<MaximalSpace>::iterator e = espacios.begin() + 1; e < espacios.end(); ++e) {
                        GenerarOpciones_Multidrop_Presion(opciones, *e);
                        if (opciones.size() > 0) {
                            vector<Bloque>::iterator mejorBloque = opciones.begin();
                            SeleccionarBloque_MaxVol(opciones, mejorBloque);
                            Empacar_Multidrop_Presion(mejorBloque, e);
                            terminar = false;
                            break;
                        }
                    }
                    if (terminar) break;
                }
            }
        }
    }
    void ConstructivoDeterministico_Multidrop_Presion() {
        while ((int)espacios.size() > 0 && hayCajasPorEmpacar) { // Falta poner restricción si hay algún tipo con caja

            // Se selecciona el espacio maximal

            vector<MaximalSpace>::iterator mejorEspacio = espacios.begin();
            double myRand1 = dist01(generator);
            if (myRand1 < 0.33) SeleccionarEspacio_Fondo_Altura(mejorEspacio);
            else if (myRand1 < 0.67) SeleccionarEspacio_Fondo_Bajo(mejorEspacio);
            else SeleccionarEspacio_Altura_Fondo(mejorEspacio);

            // Se selecciona la pieza a empacar

            vector<Bloque> opciones;
            GenerarOpciones_Multidrop_Presion(opciones, *mejorEspacio);
            if (opciones.size() > 0) {
                vector<Bloque>::iterator mejorBloque = opciones.begin();
                SeleccionarBloque_BestFit_Deterministico(opciones, mejorBloque);
                Empacar_Multidrop_Presion(mejorBloque, mejorEspacio);
            }
            else {
                if (myRand1 < 0.33) sort(espacios.begin(), espacios.end(), Orden_Fondo_Altura);
                else if (myRand1 < 0.67) sort(espacios.begin(), espacios.end(), Orden_Fondo_Bajo);
                else sort(espacios.begin(), espacios.end(), Orden_Altura_Fondo);

                // Recorrer los demás espacios en orden

                bool terminar = true;
                for (vector<MaximalSpace>::iterator e = espacios.begin() + 1; e < espacios.end(); ++e) {
                    GenerarOpciones_Multidrop_Presion(opciones, *e);
                    if (opciones.size() > 0) {
                        vector<Bloque>::iterator mejorBloque = opciones.begin();
                        SeleccionarBloque_BestFit_Deterministico(opciones, mejorBloque);
                        Empacar_Multidrop_Presion(mejorBloque, e);
                        terminar = false;
                        break;
                    }
                }
                if (terminar) break;
            }
        }
    }

    // Quitar piezas

    void Remover0() {

        // Quitar piezas empacadas

        empacados.resize((int)ceil((double)empacados.size() * 0.7f));

        // Reestablecer cantidades de las piezas y los minD

        for (vector<BoxType>::iterator bt = boxest.begin(); bt < boxest.end(); ++bt) {
            (*bt).cantidad = (*bt).cantidad0;
            if (minDx > (*bt).minDx) minDx = (*bt).minDx;
            if (minDy > (*bt).minDy) minDy = (*bt).minDy;
            if (minDz > (*bt).minDz) minDz = (*bt).minDz;
        }

        // Encontrar nuevamete los espacios maximales

        espacios = vector<MaximalSpace>(1, MaximalSpace());
        utilizacion = 0.0;
        for (vector<PackedBox>::iterator p = empacados.begin(); p < empacados.end(); ++p) {
            --boxest[(*p).id].cantidad;
            utilizacion += (*p).volumen;
            JuntarEspacios(ActualizarEspacios(*p));
            //EliminarEspaciosMinD();
        }

        // Actualizar el valor de Dmin

        ActualizarMinD();
        EliminarEspaciosMinD();
    }
    void DeterminarClienteActual() {
        if ((int)empacados.size() == 0) {
            clienteActual = 0;
            return;
        }
        clienteActual = empacados.back().cliente;

        // Se determina si hay cajas de este cliente por empacar

        for (vector<BoxType>::reverse_iterator mybt = boxest.rbegin(); mybt != boxest.rend(); ++mybt) {
            if ((*mybt).cliente != clienteActual) break;
            else if ((*mybt).cantidad > 0) return;
        }
        for (vector<BoxType>::iterator mybt = boxest.begin() + empacados.back().id + 1; mybt < boxest.end(); ++mybt) {
            if ((*mybt).cliente != clienteActual) break;
            else if ((*mybt).cantidad > 0) return;
        }
        ++clienteActual;
    }
    void ActualizarEspaciosCapaReconstruccion() {
        int clienteBorrar = clienteActual - 1;
        int x0 = 0;
        for (vector<PackedBox>::reverse_iterator p = empacados.rbegin(); p != empacados.rend(); ++p) {
            if ((*p).cliente < clienteBorrar) break;
            else if ((*p).cliente == clienteBorrar && (*p).x2 > x0) x0 = (*p).x2;
        }
        ActualizarEspaciosCapas(x0);
    }
    void ActualizarEspaciosVisibilidadReconstruccion() {
        int clienteBorrar = clienteActual - 1;
        for (int i = (int)espacios.size() - 1; i >= 0; --i) {
            MaximalSpace e = espacios[i];
            if (e.x2 < ContenedorDimx) {
                for (vector<PackedBox>::reverse_iterator p_it = empacados.rbegin(); p_it != empacados.rend(); ++p_it) {
                    if (e.x1 < (*p_it).x2 && (*p_it).cliente <= clienteBorrar) {
                        if (!(e.z1 >= (*p_it).z2 || (*p_it).z1 >= e.z2)) {
                            if (!(e.y1 >= (*p_it).y2 || (*p_it).y1 >= e.y2)) {
                                espacios.erase(espacios.begin() + i);
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    void DeterminarSoportes() {
        if ((int)empacados.size() == 0) return;

        // Actualizar espacios maximales actuales

        for (vector<MaximalSpace>::iterator e = espacios.begin(); e < espacios.end(); ++e) {
            DeterminarIndicesCajasEspacio(e);
        }

        // Se inicializan las cajas empacadas

        vector<PackedBox> todas;
        todas.reserve((int)empacados.size());
        for (vector<PackedBox>::iterator c = empacados.begin(); c < empacados.end(); ++c) {
            (*c).centroMasax = (double)((*c).x1 + (*c).x2) / 2.0;
            (*c).centroMasay = (double)((*c).y1 + (*c).y2) / 2.0;
            (*c).masaTotal = (*c).masa;
            if ((*c).soportes.size() > 0) {
                for (vector<Soporte>::iterator s = (*c).soportes.begin(); s < (*c).soportes.end(); ++s) {
                    (*s).deltaMasa = 0;
                }
            }
            todas.push_back(*c);
        }

        // Actualizar los soportes

        sort(todas.begin(), todas.end(), [](PackedBox const& p1, PackedBox const& p2)->bool {return p1.z1 > p2.z1; });
        int ini = 0;
        int indUltimo = 0;
        int miz1 = todas.front().z1;
        while (todas[indUltimo].z1 > 0) {
            for (vector<Soporte>::iterator s = todas[indUltimo].soportes.begin(); s < todas[indUltimo].soportes.end(); ++s) {
                empacados[(*s).indPB].ActualizarCentroDeMasa(*s);
            }
            ++indUltimo;
            int zSiguiente = todas[indUltimo].z1;
            if (zSiguiente < miz1) {
                for (vector<PackedBox>::iterator mip = todas.begin() + ini; mip < todas.begin() + indUltimo; ++mip) {
                    (*mip).centroMasaEnPoligono();
                    (*mip).ActualizarSoportes();
                }
                ini = indUltimo;
                miz1 = zSiguiente;
            }
        }
    }
    void RemoverPiezas() {
        Remover0();
    }
    void RemoverPiezas_Multidrop() {
        Remover0();
        DeterminarClienteActual();
        if (r_multidrop == 1) ActualizarEspaciosVisibilidadReconstruccion();
        else if (r_multidrop == 2) ActualizarEspaciosCapaReconstruccion();
    }
    void RemoverPiezas_Presion() {
        Remover0();
        DeterminarSoportes();
    }
    void RemoverPiezas_Multidrop_Presion() {
        Remover0();
        DeterminarClienteActual();
        if (r_multidrop == 1) ActualizarEspaciosVisibilidadReconstruccion();
        else if (r_multidrop == 2) ActualizarEspaciosCapaReconstruccion();
        DeterminarSoportes();
    }

    // Operadores

    const bool operator<(Container const& otro) {
        return utilizacion < otro.utilizacion;
    }
};

// Variables globales 2

Container C0;
Container incumbente;

// Funciones

const bool ValidacionesPrevias() {

    // Inconsistencia de los parámetros globales

    if (r_maxPresionItems && r_estabilidad == 0) {
        cout << "Inconsistenacia: no se puede habilitar la restricción de presión sobre los ítems si no se tiene en cuenta la estabilidad." << endl;
        return true;
    }

    // Impresión de parámetros

    cout << "Los parámetros que se usan son los siguientes" << endl;
    if (r_multidrop == 1) cout << "Con MultiDrop - Visibilidad" << endl;
    else if (r_multidrop == 2) cout << "Con Multidrop - Capas" << endl;
    else cout << "Sin MultiDrop" << endl;
    if (r_maxPresionItems) cout << "Con presión items" << endl;
    else cout << "Sin presión items" << endl;
    if (r_estabilidad == 0) cout << "Sin soporte" << endl;
    else if (r_estabilidad == 1) cout << "Con soporte parcial" << endl;
    else cout << "Con soporte completo" << endl;
    if (r_juntarEspacios == 0) cout << "Sin juntar ni expandir espacios maximales" << endl;
    else if (r_juntarEspacios == 1) cout << "Con juntar pero sin expandir espacios maximales" << endl;
    else cout << "Con juntar y expandir espacios maximales" << endl;

    return false;
}
void EliminarCajasDimensiones(vector<BoxType>& boxest) {
    for (int i = (int)boxest.size() - 1; i >= 0; --i) {
        BoxType bt = boxest[i];
        for (int j = (int)bt.boxes.size() - 1; j >= 0; --j) {
            Box b = bt.boxes[j];
            if (b.x > ContenedorDimx || b.y > ContenedorDimy || b.z > ContenedorDimz) bt.boxes.erase(bt.boxes.begin() + j);
        }
        if ((int)bt.boxes.size() == 0) {
            boxest.erase(boxest.begin() + i);
            cout << "CUIDADO: Se eliminó completamente una pieza porque no alcanza en el contenedor" << endl;
        }
    }
}
void ReadData(string const& ins) {
    string filePath = "Instances/" + ins + ".txt";
    ifstream reader(filePath);
    vector<BoxType> bts;
    char delimitador = ' ';

    // Leer contenedor

    string line, val;
    getline(reader, line, '\n');
    stringstream iss, convertor;
    iss = stringstream(line);
    getline(iss, val, delimitador); convertor = stringstream(val);
    convertor >> ContenedorDimx;
    getline(iss, val, delimitador); convertor = stringstream(val);
    convertor >> ContenedorDimy;
    getline(iss, val, delimitador); convertor = stringstream(val);
    convertor >> ContenedorDimz;
    ContenedorVol = ContenedorDimx * ContenedorDimy * ContenedorDimz;
    getline(reader, line, '\n'); // línea que no se lee

    // Leer cajas

    if (r_multidrop > 0) {
        if (r_maxPresionItems) {
            int idCaja = 0;
            while (getline(reader, line, '\n')) {
                iss = stringstream(line);
                getline(iss, val, delimitador); convertor = stringstream(val);
                int x; convertor >> x;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rx; convertor >> rx; if (rx == 1) rx = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int y; convertor >> y;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int ry; convertor >> ry; if (ry == 1) ry = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int z; convertor >> z;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rz; convertor >> rz; if (rz == 1) rz = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int q; convertor >> q;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double m; convertor >> m;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double px; convertor >> px;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double py; convertor >> py;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double pz; convertor >> pz;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int c; convertor >> c;

                // Agregar Caja

                bts.push_back(BoxType(idCaja, q, x, y, z, rx, ry, rz, m, px, py, pz, c));
                ++idCaja;
            }
        }
        else {
            int idCaja = 0;
            while (getline(reader, line, '\n')) {
                iss = stringstream(line);
                getline(iss, val, delimitador); convertor = stringstream(val);
                int x; convertor >> x;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rx; convertor >> rx; if (rx == 1) rx = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int y; convertor >> y;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int ry; convertor >> ry; if (ry == 1) ry = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int z; convertor >> z;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rz; convertor >> rz; if (rz == 1) rz = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int q; convertor >> q;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int c; convertor >> c;

                // Agregar Caja

                bts.push_back(BoxType(idCaja, q, x, y, z, rx, ry, rz, 0, c));
                ++idCaja;
            }
        }
        EliminarCajasDimensiones(bts);
        sort(bts.begin(), bts.end(), [](BoxType const& bt1, BoxType const& bt2)->bool {return bt1.cliente < bt2.cliente; });

        // Se pone el id del cliente correcto y del tipo de caja

        int clienteReal = bts.front().cliente;
        int clienteModificado = 0;
        int idV = 0;
        for (vector<BoxType>::iterator bt = bts.begin(); bt < bts.end(); ++bt, ++idV) {
            if ((*bt).cliente != clienteReal) {
                ++clienteModificado;
                clienteReal = (*bt).cliente;
            }
            (*bt).cliente = clienteModificado;
            (*bt).ActualizarIndice(idV);
        }
    }
    else {
        if (r_maxPresionItems) {
            int idCaja = 0;
            while (getline(reader, line, '\n')) {
                iss = stringstream(line);
                getline(iss, val, delimitador); convertor = stringstream(val);
                int x; convertor >> x;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rx; convertor >> rx; if (rx == 1) rx = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int y; convertor >> y;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int ry; convertor >> ry; if (ry == 1) ry = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int z; convertor >> z;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rz; convertor >> rz; if (rz == 1) rz = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int q; convertor >> q;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double m; convertor >> m;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double px; convertor >> px;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double py; convertor >> py;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double pz; convertor >> pz;

                // Agregar Caja

                bts.push_back(BoxType(idCaja, q, x, y, z, rx, ry, rz, m, px, py, pz));
                ++idCaja;
            }
        }
        else {
            int idCaja = 0;
            while (getline(reader, line, '\n')) {
                iss = stringstream(line);
                getline(iss, val, delimitador); convertor = stringstream(val);
                int x; convertor >> x;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rx; convertor >> rx; if (rx == 1) rx = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int y; convertor >> y;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int ry; convertor >> ry; if (ry == 1) ry = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int z; convertor >> z;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rz; convertor >> rz; if (rz == 1) rz = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int q; convertor >> q;

                // Agregar Caja

                bts.push_back(BoxType(idCaja, q, x, y, z, rx, ry, rz));
                ++idCaja;
            }
        }
        EliminarCajasDimensiones(bts);

        // Actualizar id del tipo de caja

        int idV = 0;
        for (vector<BoxType>::iterator bt = bts.begin(); bt < bts.end(); ++bt, ++idV) {
            (*bt).ActualizarIndice(idV);
        }
    }
    reader.close();

    // Crear contenedor

    C0 = Container(bts, 0, 0);
}
void ReadData2(string const& ins)
{
    string filePath = "Instances/" + ins + ".txt";
    ifstream reader(filePath);
    vector<BoxType> bts;
    char delimitador = '\t';

    // Leer contenedor

    string line, val;
    getline(reader, line, '\n');
    getline(reader, line, '\n');
    stringstream iss, convertor;
    iss = stringstream(line);
    getline(iss, val, delimitador); convertor = stringstream(val);
    convertor >> ContenedorDimx;
    getline(iss, val, delimitador); convertor = stringstream(val);
    convertor >> ContenedorDimy;
    getline(iss, val, delimitador); convertor = stringstream(val);
    convertor >> ContenedorDimz;
    ContenedorVol = ContenedorDimx * ContenedorDimy * ContenedorDimz;

    // Leer cajas

    if (r_multidrop > 0) {
        if (r_maxPresionItems) {
            int idCaja = 0;
            while (getline(reader, line, '\n')) {
                iss = stringstream(line);
                getline(iss, val, delimitador);
                getline(iss, val, delimitador); convertor = stringstream(val);
                int x; convertor >> x;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rx; convertor >> rx; if (rx == 1) rx = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int y; convertor >> y;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int ry; convertor >> ry; if (ry == 1) ry = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int z; convertor >> z;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rz; convertor >> rz; if (rz == 1) rz = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int q; convertor >> q;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double m; convertor >> m;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double px; convertor >> px;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double py; convertor >> py;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double pz; convertor >> pz;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int c; convertor >> c;
                /*
                int c;
                while (getline(iss, val, delimitador)) {
                    convertor = stringstream(val);
                    convertor >> c;
                }
                */

                // Agregar Caja

                bts.push_back(BoxType(idCaja, q, x, y, z, rx, ry, rz, m, px, py, pz, c));
                ++idCaja;
            }
        }
        else {
            int idCaja = 0;
            while (getline(reader, line, '\n')) {
                iss = stringstream(line);
                getline(iss, val, delimitador);
                getline(iss, val, delimitador); convertor = stringstream(val);
                int x; convertor >> x;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rx; convertor >> rx; if (rx == 1) rx = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int y; convertor >> y;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int ry; convertor >> ry; if (ry == 1) ry = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int z; convertor >> z;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rz; convertor >> rz; if (rz == 1) rz = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int q; convertor >> q;
                int c;
                while (getline(iss, val, delimitador)) {
                    convertor = stringstream(val);
                    convertor >> c;
                }

                // Agregar Caja

                bts.push_back(BoxType(idCaja, q, x, y, z, rx, ry, rz, 0, c));
                ++idCaja;
            }
        }
        EliminarCajasDimensiones(bts);
        sort(bts.begin(), bts.end(), [](BoxType const& bt1, BoxType const& bt2)->bool {return bt1.cliente < bt2.cliente; });

        // Actualizar id del cliente y del tipo de caja

        int clienteReal = bts.front().cliente;
        int clienteModificado = 0;
        int idV = 0;
        for (vector<BoxType>::iterator bt = bts.begin(); bt < bts.end(); ++bt, ++idV) {
            if ((*bt).cliente != clienteReal) {
                ++clienteModificado;
                clienteReal = (*bt).cliente;
            }
            (*bt).cliente = clienteModificado;
            (*bt).ActualizarIndice(idV);
        }
    }
    else {
        if (r_maxPresionItems) {
            int idCaja = 0;
            while (getline(reader, line, '\n')) {
                iss = stringstream(line);
                getline(iss, val, delimitador);
                getline(iss, val, delimitador); convertor = stringstream(val);
                int x; convertor >> x;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rx; convertor >> rx; if (rx == 1) rx = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int y; convertor >> y;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int ry; convertor >> ry; if (ry == 1) ry = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int z; convertor >> z;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rz; convertor >> rz; if (rz == 1) rz = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int q; convertor >> q;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double m; convertor >> m;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double px; convertor >> px;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double py; convertor >> py;
                getline(iss, val, delimitador); convertor = stringstream(val);
                double pz; convertor >> pz;

                // Agregar Caja

                bts.push_back(BoxType(idCaja, q, x, y, z, rx, ry, rz, m, px, py, pz));
                ++idCaja;
            }
        }
        else {
            int idCaja = 0;
            while (getline(reader, line, '\n')) {
                iss = stringstream(line);
                getline(iss, val, delimitador);
                getline(iss, val, delimitador); convertor = stringstream(val);
                int x; convertor >> x;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rx; convertor >> rx; if (rx == 1) rx = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int y; convertor >> y;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int ry; convertor >> ry; if (ry == 1) ry = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int z; convertor >> z;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int rz; convertor >> rz; if (rz == 1) rz = 3;
                getline(iss, val, delimitador); convertor = stringstream(val);
                int q; convertor >> q;

                // Agregar Caja

                bts.push_back(BoxType(idCaja, q, x, y, z, rx, ry, rz));
                ++idCaja;
            }
        }
        EliminarCajasDimensiones(bts);

        // Actualizar id del tipo de caja

        int idV = 0;
        for (vector<BoxType>::iterator bt = bts.begin(); bt < bts.end(); ++bt, ++idV) {
            (*bt).ActualizarIndice(idV);
        }
    }
    reader.close();

    // Crear contenedor

    C0 = Container(bts, 0, 0);
}
const bool ValidacionesDatos() {
    if (ContenedorDimx <= 0 || ContenedorDimy <= 0 || ContenedorDimz <= 0) {
        cout << "Alguna dimension del contenedor es menor o igual a cero" << endl;
        return true;
    }
    for (vector<BoxType>::iterator bt = C0.boxest.begin(); bt < C0.boxest.end(); ++bt) {
        if ((*bt).cantidad < 0) {
            cout << "Alguna cantidad es menor a cero" << endl;
            return true;
        }
        if ((*bt).masa < 0) {
            cout << "Alguna masa es menor a cero" << endl;
            return true;
        }
        Box b = (*bt).boxes.front();
        if (b.x <= 0 || b.y <= 0 || b.z <= 0) {
            cout << "Alguna caja tiene una dimensión menor o igual a cero" << endl;
            return true;
        }
        for (vector<Box>::iterator bb = (*bt).boxes.begin(); bb < (*bt).boxes.end(); ++bb) {
            if ((*bb).presion < 0) {
                cout << "Alguna presión es menor a cero" << endl;
                return true;
            }
        }
    }
    return false;
}
void UpdateTime()
{
    auto miDuracion = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tIni);
    duracion = (double)miDuracion.count() / 1000.0;
}
void UpdateTimeParal(double& paralDur)
{
    auto miDuracion = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - tIni);
    paralDur = (double)miDuracion.count() / 1000.0;
}
void WriteData(string const& ins, int const& repeticion, int const& totalIteraciones, int th, int dur)
{
    // Tiempo

    UpdateTime();
    /*
    for (vector<BoxType>::iterator b_it = incumbente.boxest.begin(); b_it < incumbente.boxest.end(); ++b_it) {
        if ((*b_it).cantidad > 0) {
            cout << "Id = " << (*b_it).id << " cliente = " << (*b_it).cliente << " Cantidad = " << (*b_it).cantidad << " idReal = " << (*b_it).idReal << endl;
        }
    }
    */
    // Escritura del archivo

    string filePath = "Solutions/R_" + ins;
    filePath += "_MD" + to_string(r_multidrop);
    if (r_maxPresionItems) filePath += "_P1";
    else filePath += "_P0";
    filePath += "_S" + to_string(r_estabilidad);
    filePath += "_EM" + to_string(r_juntarEspacios);
    filePath += "_st_rep" + to_string(repeticion);
    filePath += "_nt" + to_string(th);
    filePath += "_d" + to_string(dur);
    filePath += ".txt";
    ofstream writer(filePath);
    writer << "Utilizacion " << to_string(incumbente.utilizacion) << " VolCajasEmpacado " << to_string(incumbente.volEmpacado) << " Duracion[s] " << to_string(duracion) << " Iteraciones " << to_string(totalIteraciones) << endl;
    int id = 0;
    for (vector<PackedBox>::iterator p = incumbente.empacados.begin(); p < incumbente.empacados.end(); ++p, ++id) {
        writer << to_string(id) << " " << to_string((*p).id) << " " << to_string((*p).cliente) << " " + to_string((*p).x1) << " " << to_string((*p).y1) << " " << to_string((*p).z1) << " " << to_string((*p).x2) << " " << to_string((*p).y2) << " " << to_string((*p).z2) << endl;
    }
    writer.close();
    cout << "ok" << endl;
}
vector<Alpha>::iterator DeterminarAlpha() {
    sort(alphas.begin(), alphas.end());
    double val = dist01(generator);
    vector<double>::iterator p_it = pesos_acum.begin();
    for (vector<Alpha>::iterator a = alphas.begin(); a < alphas.end(); ++a, ++p_it) {
        if (val < *p_it) {
            return a;
        }
    }
    return alphas.end() - 1;
}
vector<Alpha>::iterator DeterminarAlphaParal(vector<Alpha>& vecAlphas) {
    sort(vecAlphas.begin(), vecAlphas.end());
    double val = dist01(generator);
    vector<double>::iterator p_it = pesos_acum.begin();
    for (vector<Alpha>::iterator a = vecAlphas.begin(); a < vecAlphas.end(); ++a, ++p_it) {
        if (val < *p_it) {
            return a;
        }
    }
    return vecAlphas.end() - 1;
}
void ActualizarAlpha(Alpha& a, double const& utilizacion) {
    a.Actualizar(utilizacion);
}
int ActualizarIncumbente(Container& C) {
    if (C.utilizacion > incumbente.utilizacion) {
        incumbente = Container(C, false);
        if (!incumbente.hayCajasPorEmpacar) return -1;
        return 1;
    }
    return 0;
}
void ReiniciarVariables() {
    alphas.clear();
    alphas = vector<Alpha>({ Alpha(0.1), Alpha(0.2), Alpha(0.3), Alpha(0.4), Alpha(0.5), Alpha(0.6), Alpha(0.7), Alpha(0.8), Alpha(0.9), Alpha(1.0) });
    tIni = chrono::high_resolution_clock::now();
    duracion = 0;
    //generator.seed(chrono::system_clock::now().time_since_epoch().count());
    generator.seed(43);
}

// Main

int main(int argc, char** argv) {

    // Parámetros

    string ins = "";
    int maxTime = 60;
    int nThreads = 1;
    for (int i = 1; i < argc - 1; i += 2) {
        if (argc - 1 >= i + 1) {
            if (string(argv[i]) == "-ins") ins = argv[i + 1];
            else if (string(argv[i]) == "-nThreads") nThreads = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-maxTime") maxTime = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-estabilidad") r_estabilidad = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-juntar") r_juntarEspacios = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-presion") {
                if (string(argv[i + 1]) == "1") r_maxPresionItems = true;
                else r_maxPresionItems = false;
            }
            else if (string(argv[i]) == "-multidrop") r_multidrop = atoi(argv[i + 1]);
            else {
                cout << "Mal en par�metros" << endl;
                return 0;
            }
        }
    }
    omp_set_dynamic(0);
    omp_set_num_threads(nThreads);

    // Ciclo de repeticiones

    ReiniciarVariables();

    // Leer datos

    ReadData2(ins);
    incumbente = Container(C0, false);
    errorArea = 0.5 / (double)(ContenedorDimx * ContenedorDimy);
    double totalPesoAcum = 0;
    int j = 1;
    for (vector<double>::iterator p_it = pesos_acum.begin(); p_it < pesos_acum.end(); ++p_it, ++j) {
        (*p_it) = 1.0 / (double)j;
        totalPesoAcum += *p_it;
    }
    for (vector<double>::iterator p_it = pesos_acum.begin(); p_it < pesos_acum.end(); ++p_it) {
        (*p_it) /= totalPesoAcum;
    }
    for (vector<double>::iterator p_it = pesos_acum.begin() + 1; p_it < pesos_acum.end(); ++p_it) {
        (*p_it) += *(p_it - 1);
    }

    // Volumen total

    totalVolumen = 0;
    for (vector<BoxType>::iterator b_it = C0.boxest.begin(); b_it < C0.boxest.end(); ++b_it) {
        totalVolumen += (*b_it).cantidad0 * (*b_it).volumenReal;
    }
    double MaxUtilizacionPosible = totalVolumen / ContenedorVol * 100.0;
    totalVolumen /= 100.0;
    if (r_multidrop > 0 && totalVolumen > ContenedorVol) { // Se eliminan los clientes que no alcanzan a empacarse por el volumen de los anteriores
        int miVol = 0;
        vector<BoxType>::iterator b_it = C0.boxest.begin();
        int clienteHasta = 0;
        for (; b_it < C0.boxest.end(); ++b_it) {
            miVol += (*b_it).cantidad0 * (*b_it).volumenReal;
            if (miVol >= ContenedorVol) {
                clienteHasta = (*b_it).cliente;
                break;
            }
        }
        for (; b_it < C0.boxest.end(); ++b_it) {
            if ((*b_it).cliente > clienteHasta) {
                C0.boxest.resize(distance(C0.boxest.begin(), b_it - 1));
                break;
            }
        }
    }

    // Variables de paralelización

    vector<int> paralNIter(nThreads, 0);
    vector<int> paralThreshold_iter(nThreads, 0);
    list<int> paralThreshold_val({ 0 });
    vector<int> paralHolguraMultiplicador(nThreads, 1);
    vector<vector<Alpha>> paralAlphas(nThreads, vector<Alpha>({ Alpha(0.1), Alpha(0.2), Alpha(0.3), Alpha(0.4), Alpha(0.5), Alpha(0.6), Alpha(0.7), Alpha(0.8), Alpha(0.9), Alpha(1.0) }));
    vector<Container> paralIncumbente(nThreads, C0);
    vector<double> paralDuraciones(nThreads, 0);

    // Adquisición de datos iniciales

    if (r_multidrop > 0) {
        if (r_maxPresionItems) {
#pragma omp parallel for
            for (int i = 0; i < alphas.size(); ++i) {
                int idThread = omp_get_thread_num();
                if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                    ++paralNIter[idThread];
                    Alpha& a = paralAlphas[idThread][i];
                    Container C(C0, true);
                    C.alpha = a.val;
                    C.Constructivo_Multidrop_Presion();
                    ActualizarAlpha(a, C.utilizacion);
                    if (holi == 16) {
                        incumbente = C;
                        WriteData(ins + "PRUEBA", 0, (int)accumulate(paralNIter.begin(), paralNIter.end(), 0), nThreads, (int)maxTime);
                        return 1;
                    }
                    int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                    if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                        paralIncumbente[idThread] = Container(C, false);
                        paralThreshold_val.push_back(C.utilizacion);
                    }
                    if (C.utilizacion > BestUtil) {

                        // Búsqueda local

                        C.RemoverPiezas_Multidrop_Presion();
                        C.ConstructivoDeterministico_Multidrop_Presion();
                        ActualizarAlpha(a, C.utilizacion);
                        if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                            paralIncumbente[idThread] = Container(C, false);
                            paralThreshold_val.push_back(C.utilizacion);
                        }
                        paralThreshold_iter[idThread] = 0;
                    }
                    else {
                        ++paralThreshold_iter[idThread];
                        if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                            paralThreshold_iter[idThread] = 0;
                            paralHolguraMultiplicador[idThread] *= 0.8;
                        }
                    }

                    // Parar por tiempo

                    UpdateTimeParal(paralDuraciones[idThread]);
                }
            }
        }
        else {
#pragma omp parallel for
            for (int i = 0; i < alphas.size(); ++i) {
                int idThread = omp_get_thread_num();
                if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                    ++paralNIter[idThread];
                    Alpha& a = paralAlphas[idThread][i];
                    Container C(C0, true);
                    C.alpha = a.val;
                    C.Constructivo_Multidrop();
                    ActualizarAlpha(a, C.utilizacion);
                    int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                    if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                        paralIncumbente[idThread] = Container(C, false);
                        paralThreshold_val.push_back(C.utilizacion);
                    }
                    if (C.utilizacion > BestUtil) {

                        // Búsqueda local

                        C.RemoverPiezas_Multidrop();
                        C.ConstructivoDeterministico_Multidrop();
                        ActualizarAlpha(a, C.utilizacion);
                        if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                            paralIncumbente[idThread] = Container(C, false);
                            paralThreshold_val.push_back(C.utilizacion);
                        }
                        paralThreshold_iter[idThread] = 0;
                    }
                    else {
                        ++paralThreshold_iter[idThread];
                        if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                            paralThreshold_iter[idThread] = 0;
                            paralHolguraMultiplicador[idThread] *= 0.8;
                        }
                    }

                    // Parar por tiempo

                    if (idThread == 0) {
                        UpdateTimeParal(paralDuraciones[idThread]);
                    }
                }
            }
        }
    }
    else {
        if (r_maxPresionItems) {
#pragma omp parallel for
            for (int i = 0; i < alphas.size(); ++i) {
                int idThread = omp_get_thread_num();
                if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                    ++paralNIter[idThread];
                    Alpha& a = paralAlphas[idThread][i];
                    Container C(C0, true);
                    C.alpha = a.val;
                    C.Constructivo_Presion();
                    ActualizarAlpha(a, C.utilizacion);
                    int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                    if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                        paralIncumbente[idThread] = Container(C, false);
                        paralThreshold_val.push_back(C.utilizacion);
                    }
                    if (C.utilizacion > BestUtil) {

                        // Búsqueda local

                        C.RemoverPiezas_Presion();
                        C.ConstructivoDeterministico_Presion();
                        ActualizarAlpha(a, C.utilizacion);
                        if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                            paralIncumbente[idThread] = Container(C, false);
                            paralThreshold_val.push_back(C.utilizacion);
                        }
                        paralThreshold_iter[idThread] = 0;
                    }
                    else {
                        ++paralThreshold_iter[idThread];
                        if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                            paralThreshold_iter[idThread] = 0;
                            paralHolguraMultiplicador[idThread] *= 0.8;
                        }
                    }

                    // Parar por tiempo

                    if (idThread == 0) {
                        UpdateTimeParal(paralDuraciones[idThread]);
                    }
                }
            }
        }
        else {
#pragma omp parallel for
            for (int i = 0; i < alphas.size(); ++i) {
                int idThread = omp_get_thread_num();
                if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                    ++paralNIter[idThread];
                    Alpha& a = paralAlphas[idThread][i];
                    Container C(C0, true);
                    C.alpha = a.val;
                    C.Constructivo();
                    ActualizarAlpha(a, C.utilizacion);
                    int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                    if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                        paralIncumbente[idThread] = Container(C, false);
                        paralThreshold_val.push_back(C.utilizacion);
                    }
                    if (C.utilizacion > BestUtil) {

                        // Búsqueda local

                        C.RemoverPiezas();
                        C.ConstructivoDeterministico();
                        ActualizarAlpha(a, C.utilizacion);
                        if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                            paralIncumbente[idThread] = Container(C, false);
                            paralThreshold_val.push_back(C.utilizacion);
                        }
                        paralThreshold_iter[idThread] = 0;
                    }
                    else {
                        ++paralThreshold_iter[idThread];
                        if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                            paralThreshold_iter[idThread] = 0;
                            paralHolguraMultiplicador[idThread] *= 0.8;
                        }
                    }

                    // Parar por tiempo

                    if (idThread == 0) {
                        UpdateTimeParal(paralDuraciones[idThread]);
                    }
                }
            }
        }
    }

    // Resto de la metaheurística

    UpdateTime();
    incumbente = *max_element(paralIncumbente.begin(), paralIncumbente.end());
    if (duracion >= maxTime || !incumbente.hayCajasPorEmpacar) {
        WriteData(ins, 0, (int)accumulate(paralNIter.begin(), paralNIter.end(), 0), nThreads, (int)maxTime);
        return 1;
    }
    int maxIter = 100000;
    if (r_multidrop > 0) {
        if (r_maxPresionItems) {
            while (duracion < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
#pragma omp parallel for
                for (int i = 0; i < maxIter; ++i) {
                    int idThread = omp_get_thread_num();
                    if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                        ++paralNIter[idThread];
                        vector<Alpha>::iterator a = DeterminarAlphaParal(paralAlphas[idThread]);
                        Container C(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo_Multidrop_Presion();
                        ActualizarAlpha(*a, C.utilizacion);
                        int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                        if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                            paralIncumbente[idThread] = Container(C, false);
                            paralThreshold_val.push_back(C.utilizacion);
                        }
                        if (C.utilizacion > BestUtil) {

                            // Búsqueda local

                            C.RemoverPiezas_Multidrop_Presion();
                            C.ConstructivoDeterministico_Multidrop_Presion();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                paralIncumbente[idThread] = Container(C, false);
                                paralThreshold_val.push_back(C.utilizacion);
                            }
                            paralThreshold_iter[idThread] = 0;
                        }
                        else {
                            ++paralThreshold_iter[idThread];
                            if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                paralThreshold_iter[idThread] = 0;
                                paralHolguraMultiplicador[idThread] *= 0.8;
                            }
                        }

                        // Parar por tiempo

                        UpdateTimeParal(paralDuraciones[idThread]);
                    }
                }
                UpdateTime();
            }
        }
        else {
            while (duracion < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
#pragma omp parallel for
                for (int i = 0; i < maxIter; ++i) {
                    int idThread = omp_get_thread_num();
                    if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                        ++paralNIter[idThread];
                        vector<Alpha>::iterator a = DeterminarAlphaParal(paralAlphas[idThread]);
                        Container C(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo_Multidrop();
                        ActualizarAlpha(*a, C.utilizacion);
                        int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                        if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                            paralIncumbente[idThread] = Container(C, false);
                            paralThreshold_val.push_back(C.utilizacion);
                        }
                        if (C.utilizacion > BestUtil) {

                            // Búsqueda local

                            C.RemoverPiezas_Multidrop();
                            C.ConstructivoDeterministico_Multidrop();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                paralIncumbente[idThread] = Container(C, false);
                                paralThreshold_val.push_back(C.utilizacion);
                            }
                            paralThreshold_iter[idThread] = 0;
                        }
                        else {
                            ++paralThreshold_iter[idThread];
                            if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                paralThreshold_iter[idThread] = 0;
                                paralHolguraMultiplicador[idThread] *= 0.8;
                            }
                        }

                        // Parar por tiempo

                        UpdateTimeParal(paralDuraciones[idThread]);
                    }
                }
                UpdateTime();
            }
        }
    }
    else {
        if (r_maxPresionItems) {
            while (duracion < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
#pragma omp parallel for
                for (int i = 0; i < maxIter; ++i) {
                    int idThread = omp_get_thread_num();
                    if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                        ++paralNIter[idThread];
                        vector<Alpha>::iterator a = DeterminarAlphaParal(paralAlphas[idThread]);
                        Container C(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo_Presion();
                        ActualizarAlpha(*a, C.utilizacion);
                        int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                        if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                            paralIncumbente[idThread] = Container(C, false);
                            paralThreshold_val.push_back(C.utilizacion);
                        }
                        if (C.utilizacion > BestUtil) {

                            // Búsqueda local

                            C.RemoverPiezas_Presion();
                            C.ConstructivoDeterministico_Presion();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                paralIncumbente[idThread] = Container(C, false);
                                paralThreshold_val.push_back(C.utilizacion);
                            }
                            paralThreshold_iter[idThread] = 0;
                        }
                        else {
                            ++paralThreshold_iter[idThread];
                            if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                paralThreshold_iter[idThread] = 0;
                                paralHolguraMultiplicador[idThread] *= 0.8;
                            }
                        }

                        // Parar por tiempo

                        UpdateTimeParal(paralDuraciones[idThread]);
                    }
                }
                UpdateTime();
            }
        }
        else {
            while (duracion < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
#pragma omp parallel for
                for (int i = 0; i < maxIter; ++i) {
                    int idThread = omp_get_thread_num();
                    if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                        ++paralNIter[idThread];
                        vector<Alpha>::iterator a = DeterminarAlphaParal(paralAlphas[idThread]);
                        Container C(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo();
                        ActualizarAlpha(*a, C.utilizacion);
                        int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                        if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                            paralIncumbente[idThread] = Container(C, false);
                            paralThreshold_val.push_back(C.utilizacion);
                        }
                        if (C.utilizacion > BestUtil) {

                            // Búsqueda local

                            C.RemoverPiezas();
                            C.ConstructivoDeterministico();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                paralIncumbente[idThread] = Container(C, false);
                                paralThreshold_val.push_back(C.utilizacion);
                            }
                            paralThreshold_iter[idThread] = 0;
                        }
                        else {
                            ++paralThreshold_iter[idThread];
                            if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                paralThreshold_iter[idThread] = 0;
                                paralHolguraMultiplicador[idThread] *= 0.8;
                            }
                        }

                        // Parar por tiempo

                        UpdateTimeParal(paralDuraciones[idThread]);
                    }
                }
                UpdateTime();
            }
        }
    }

    // Escribir resultados

    incumbente = *max_element(paralIncumbente.begin(), paralIncumbente.end());
    WriteData(ins, 0, (int)accumulate(paralNIter.begin(), paralNIter.end(), 0), nThreads, (int)maxTime);
    return 0;
}

// Multi Thread
/*
int main(int argc, char** argv) {

    // Parámetros

    string ins = "";
    int maxTime = 60;
    int nThreads = 1;
    int iterIni = 1;
    int iterFin = 1;// Inclusivo
    int insIni = 0;
    int insFin = 0; // Inclusivo
    for (int i = 1; i < argc - 1; i += 2) {
        if (argc - 1 >= i + 1) {
            if (string(argv[i]) == "-ins") ins = argv[i + 1];
            else if (string(argv[i]) == "-nThreads") nThreads = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-maxTime") maxTime = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-estabilidad") r_estabilidad = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-juntar") r_juntarEspacios = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-ini") iterIni = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-fin") iterFin = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-insIni") insIni = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-insFin") insFin = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-presion") {
                if (string(argv[i + 1]) == "1") r_maxPresionItems = true;
                else r_maxPresionItems = false;
            }
            else if (string(argv[i]) == "-multidrop") r_multidrop = atoi(argv[i + 1]);
            else {
                cout << "Mal en par�metros" << endl;
                return 0;
            }
        }
    }
    omp_set_dynamic(0);
    omp_set_num_threads(nThreads);

    // Determinar los clientes

    string ins0 = ins;
    vector<string> misInstancias({ "" });
    if (insFin - insIni >= 1) {
        misInstancias.clear();
        for (int i = insIni; i <= insFin; ++i) {
            misInstancias.push_back(to_string(i));
        }
    }

    // Ciclo de repeticiones

    for (vector<string>::iterator s_it = misInstancias.begin(); s_it < misInstancias.end(); ++s_it) {
        ins = ins0;
        if ((*s_it) != "") {
            ins.replace(ins.find("__"), 2, "_" + *s_it + "_");
        }
        for (int repeticion = iterIni; repeticion <= iterFin; ++repeticion) {
            ReiniciarVariables();

            // Leer datos

            ReadData2(ins);
            incumbente = Container(C0, false);
            errorArea = 0.5 / (double)(ContenedorDimx * ContenedorDimy);
            double totalPesoAcum = 0;
            int j = 1;
            for (vector<double>::iterator p_it = pesos_acum.begin(); p_it < pesos_acum.end(); ++p_it, ++j) {
                (*p_it) = 1.0 / (double)j;
                totalPesoAcum += *p_it;
            }
            for (vector<double>::iterator p_it = pesos_acum.begin(); p_it < pesos_acum.end(); ++p_it) {
                (*p_it) /= totalPesoAcum;
            }
            for (vector<double>::iterator p_it = pesos_acum.begin() + 1; p_it < pesos_acum.end(); ++p_it) {
                (*p_it) += *(p_it - 1);
            }

            // Volumen total

            totalVolumen = 0;
            for (vector<BoxType>::iterator b_it = C0.boxest.begin(); b_it < C0.boxest.end(); ++b_it) {
                totalVolumen += (*b_it).cantidad0 * (*b_it).volumenReal;
            }
            double MaxUtilizacionPosible = totalVolumen / ContenedorVol * 100.0;
            totalVolumen /= 100.0;
            if (r_multidrop > 0 && totalVolumen > ContenedorVol) { // Se eliminan los clientes que no alcanzan a empacarse por el volumen de los anteriores
                int miVol = 0;
                vector<BoxType>::iterator b_it = C0.boxest.begin();
                int clienteHasta = 0;
                for (; b_it < C0.boxest.end(); ++b_it) {
                    miVol += (*b_it).cantidad0 * (*b_it).volumenReal;
                    if (miVol >= ContenedorVol) {
                        clienteHasta = (*b_it).cliente;
                        break;
                    }
                }
                for (; b_it < C0.boxest.end(); ++b_it) {
                    if ((*b_it).cliente > clienteHasta) {
                        C0.boxest.resize(distance(C0.boxest.begin(), b_it - 1));
                        break;
                    }
                }
            }

            // Variables de paralelización

            vector<int> paralNIter(nThreads, 0);
            vector<int> paralThreshold_iter(nThreads, 0);
            list<int> paralThreshold_val({ 0 });
            vector<int> paralHolguraMultiplicador(nThreads, 1);
            vector<vector<Alpha>> paralAlphas(nThreads, vector<Alpha>({ Alpha(0.1), Alpha(0.2), Alpha(0.3), Alpha(0.4), Alpha(0.5), Alpha(0.6), Alpha(0.7), Alpha(0.8), Alpha(0.9), Alpha(1.0) }));
            vector<Container> paralIncumbente(nThreads, C0);
            vector<double> paralDuraciones(nThreads, 0);

            // Adquisición de datos iniciales

            if (r_multidrop > 0) {
                if (r_maxPresionItems) {
#pragma omp parallel for
                    for (int i = 0; i < alphas.size(); ++i) {
                        int idThread = omp_get_thread_num();
                        if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                            ++paralNIter[idThread];
                            Alpha& a = paralAlphas[idThread][i];
                            Container C(C0, true);
                            C.alpha = a.val;
                            C.Constructivo_Multidrop_Presion();
                            ActualizarAlpha(a, C.utilizacion);
                            int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                            if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                paralIncumbente[idThread] = Container(C, false);
                                paralThreshold_val.push_back(C.utilizacion);
                            }
                            if (C.utilizacion > BestUtil) {

                                // Búsqueda local

                                C.RemoverPiezas_Multidrop_Presion();
                                C.ConstructivoDeterministico_Multidrop_Presion();
                                ActualizarAlpha(a, C.utilizacion);
                                if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                    paralIncumbente[idThread] = Container(C, false);
                                    paralThreshold_val.push_back(C.utilizacion);
                                }
                                paralThreshold_iter[idThread] = 0;
                            }
                            else {
                                ++paralThreshold_iter[idThread];
                                if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                    paralThreshold_iter[idThread] = 0;
                                    paralHolguraMultiplicador[idThread] *= 0.8;
                                }
                            }

                            // Parar por tiempo

                            UpdateTimeParal(paralDuraciones[idThread]);
                        }
                    }
                }
                else {
#pragma omp parallel for
                    for (int i = 0; i < alphas.size(); ++i) {
                        int idThread = omp_get_thread_num();
                        if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                            ++paralNIter[idThread];
                            Alpha& a = paralAlphas[idThread][i];
                            Container C(C0, true);
                            C.alpha = a.val;
                            C.Constructivo_Multidrop();
                            ActualizarAlpha(a, C.utilizacion);
                            int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                            if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                paralIncumbente[idThread] = Container(C, false);
                                paralThreshold_val.push_back(C.utilizacion);
                            }
                            if (C.utilizacion > BestUtil) {

                                // Búsqueda local

                                C.RemoverPiezas_Multidrop();
                                C.ConstructivoDeterministico_Multidrop();
                                ActualizarAlpha(a, C.utilizacion);
                                if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                    paralIncumbente[idThread] = Container(C, false);
                                    paralThreshold_val.push_back(C.utilizacion);
                                }
                                paralThreshold_iter[idThread] = 0;
                            }
                            else {
                                ++paralThreshold_iter[idThread];
                                if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                    paralThreshold_iter[idThread] = 0;
                                    paralHolguraMultiplicador[idThread] *= 0.8;
                                }
                            }

                            // Parar por tiempo

                            if (idThread == 0) {
                                UpdateTimeParal(paralDuraciones[idThread]);
                            }
                        }
                    }
                }
            }
            else {
                if (r_maxPresionItems) {
#pragma omp parallel for
                    for (int i = 0; i < alphas.size(); ++i) {
                        int idThread = omp_get_thread_num();
                        if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                            ++paralNIter[idThread];
                            Alpha& a = paralAlphas[idThread][i];
                            Container C(C0, true);
                            C.alpha = a.val;
                            C.Constructivo_Presion();
                            ActualizarAlpha(a, C.utilizacion);
                            int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                            if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                paralIncumbente[idThread] = Container(C, false);
                                paralThreshold_val.push_back(C.utilizacion);
                            }
                            if (C.utilizacion > BestUtil) {

                                // Búsqueda local

                                C.RemoverPiezas_Presion();
                                C.ConstructivoDeterministico_Presion();
                                ActualizarAlpha(a, C.utilizacion);
                                if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                    paralIncumbente[idThread] = Container(C, false);
                                    paralThreshold_val.push_back(C.utilizacion);
                                }
                                paralThreshold_iter[idThread] = 0;
                            }
                            else {
                                ++paralThreshold_iter[idThread];
                                if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                    paralThreshold_iter[idThread] = 0;
                                    paralHolguraMultiplicador[idThread] *= 0.8;
                                }
                            }

                            // Parar por tiempo

                            if (idThread == 0) {
                                UpdateTimeParal(paralDuraciones[idThread]);
                            }
                        }
                    }
                }
                else {
#pragma omp parallel for
                    for (int i = 0; i < alphas.size(); ++i) {
                        int idThread = omp_get_thread_num();
                        if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                            ++paralNIter[idThread];
                            Alpha& a = paralAlphas[idThread][i];
                            Container C(C0, true);
                            C.alpha = a.val;
                            C.Constructivo();
                            ActualizarAlpha(a, C.utilizacion);
                            int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                            if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                paralIncumbente[idThread] = Container(C, false);
                                paralThreshold_val.push_back(C.utilizacion);
                            }
                            if (C.utilizacion > BestUtil) {

                                // Búsqueda local

                                C.RemoverPiezas();
                                C.ConstructivoDeterministico();
                                ActualizarAlpha(a, C.utilizacion);
                                if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                    paralIncumbente[idThread] = Container(C, false);
                                    paralThreshold_val.push_back(C.utilizacion);
                                }
                                paralThreshold_iter[idThread] = 0;
                            }
                            else {
                                ++paralThreshold_iter[idThread];
                                if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                    paralThreshold_iter[idThread] = 0;
                                    paralHolguraMultiplicador[idThread] *= 0.8;
                                }
                            }

                            // Parar por tiempo

                            if (idThread == 0) {
                                UpdateTimeParal(paralDuraciones[idThread]);
                            }
                        }
                    }
                }
            }

            // Resto de la metaheurística

            UpdateTime();
            incumbente = *max_element(paralIncumbente.begin(), paralIncumbente.end());
            if (duracion >= maxTime || !incumbente.hayCajasPorEmpacar) {
                WriteData(ins, repeticion, (int)accumulate(paralNIter.begin(), paralNIter.end(), 0), nThreads, (int)maxTime);
                continue;
            }
            int maxIter = 100000;
            if (r_multidrop > 0) {
                if (r_maxPresionItems) {
                    while (duracion < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
#pragma omp parallel for
                        for (int i = 0; i < maxIter; ++i) {
                            int idThread = omp_get_thread_num();
                            if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                                ++paralNIter[idThread];
                                vector<Alpha>::iterator a = DeterminarAlphaParal(paralAlphas[idThread]);
                                Container C(C0, true);
                                C.alpha = (*a).val;
                                C.Constructivo_Multidrop_Presion();
                                ActualizarAlpha(*a, C.utilizacion);
                                int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                                if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                    paralIncumbente[idThread] = Container(C, false);
                                    paralThreshold_val.push_back(C.utilizacion);
                                }
                                if (C.utilizacion > BestUtil) {

                                    // Búsqueda local

                                    C.RemoverPiezas_Multidrop_Presion();
                                    C.ConstructivoDeterministico_Multidrop_Presion();
                                    ActualizarAlpha(*a, C.utilizacion);
                                    if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                        paralIncumbente[idThread] = Container(C, false);
                                        paralThreshold_val.push_back(C.utilizacion);
                                    }
                                    paralThreshold_iter[idThread] = 0;
                                }
                                else {
                                    ++paralThreshold_iter[idThread];
                                    if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                        paralThreshold_iter[idThread] = 0;
                                        paralHolguraMultiplicador[idThread] *= 0.8;
                                    }
                                }

                                // Parar por tiempo

                                UpdateTimeParal(paralDuraciones[idThread]);
                            }
                        }
                        UpdateTime();
                    }
                }
                else {
                    while (duracion < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
#pragma omp parallel for
                        for (int i = 0; i < maxIter; ++i) {
                            int idThread = omp_get_thread_num();
                            if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                                ++paralNIter[idThread];
                                vector<Alpha>::iterator a = DeterminarAlphaParal(paralAlphas[idThread]);
                                Container C(C0, true);
                                C.alpha = (*a).val;
                                C.Constructivo_Multidrop();
                                ActualizarAlpha(*a, C.utilizacion);
                                int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                                if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                    paralIncumbente[idThread] = Container(C, false);
                                    paralThreshold_val.push_back(C.utilizacion);
                                }
                                if (C.utilizacion > BestUtil) {

                                    // Búsqueda local

                                    C.RemoverPiezas_Multidrop();
                                    C.ConstructivoDeterministico_Multidrop();
                                    ActualizarAlpha(*a, C.utilizacion);
                                    if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                        paralIncumbente[idThread] = Container(C, false);
                                        paralThreshold_val.push_back(C.utilizacion);
                                    }
                                    paralThreshold_iter[idThread] = 0;
                                }
                                else {
                                    ++paralThreshold_iter[idThread];
                                    if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                        paralThreshold_iter[idThread] = 0;
                                        paralHolguraMultiplicador[idThread] *= 0.8;
                                    }
                                }

                                // Parar por tiempo

                                UpdateTimeParal(paralDuraciones[idThread]);
                            }
                        }
                        UpdateTime();
                    }
                }
            }
            else {
                if (r_maxPresionItems) {
                    while (duracion < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
#pragma omp parallel for
                        for (int i = 0; i < maxIter; ++i) {
                            int idThread = omp_get_thread_num();
                            if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                                ++paralNIter[idThread];
                                vector<Alpha>::iterator a = DeterminarAlphaParal(paralAlphas[idThread]);
                                Container C(C0, true);
                                C.alpha = (*a).val;
                                C.Constructivo_Presion();
                                ActualizarAlpha(*a, C.utilizacion);
                                int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                                if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                    paralIncumbente[idThread] = Container(C, false);
                                    paralThreshold_val.push_back(C.utilizacion);
                                }
                                if (C.utilizacion > BestUtil) {

                                    // Búsqueda local

                                    C.RemoverPiezas_Presion();
                                    C.ConstructivoDeterministico_Presion();
                                    ActualizarAlpha(*a, C.utilizacion);
                                    if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                        paralIncumbente[idThread] = Container(C, false);
                                        paralThreshold_val.push_back(C.utilizacion);
                                    }
                                    paralThreshold_iter[idThread] = 0;
                                }
                                else {
                                    ++paralThreshold_iter[idThread];
                                    if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                        paralThreshold_iter[idThread] = 0;
                                        paralHolguraMultiplicador[idThread] *= 0.8;
                                    }
                                }

                                // Parar por tiempo

                                UpdateTimeParal(paralDuraciones[idThread]);
                            }
                        }
                        UpdateTime();
                    }
                }
                else {
                    while (duracion < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
#pragma omp parallel for
                        for (int i = 0; i < maxIter; ++i) {
                            int idThread = omp_get_thread_num();
                            if (paralDuraciones[idThread] < maxTime && *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) < MaxUtilizacionPosible) {
                                ++paralNIter[idThread];
                                vector<Alpha>::iterator a = DeterminarAlphaParal(paralAlphas[idThread]);
                                Container C(C0, true);
                                C.alpha = (*a).val;
                                C.Constructivo();
                                ActualizarAlpha(*a, C.utilizacion);
                                int BestUtil = *max_element(paralThreshold_val.begin(), paralThreshold_val.end()) * paralHolguraMultiplicador[idThread];
                                if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                    paralIncumbente[idThread] = Container(C, false);
                                    paralThreshold_val.push_back(C.utilizacion);
                                }
                                if (C.utilizacion > BestUtil) {

                                    // Búsqueda local

                                    C.RemoverPiezas();
                                    C.ConstructivoDeterministico();
                                    ActualizarAlpha(*a, C.utilizacion);
                                    if (C.utilizacion > paralIncumbente[idThread].utilizacion) {
                                        paralIncumbente[idThread] = Container(C, false);
                                        paralThreshold_val.push_back(C.utilizacion);
                                    }
                                    paralThreshold_iter[idThread] = 0;
                                }
                                else {
                                    ++paralThreshold_iter[idThread];
                                    if (paralThreshold_iter[idThread] >= threshold_MaxIter) {
                                        paralThreshold_iter[idThread] = 0;
                                        paralHolguraMultiplicador[idThread] *= 0.8;
                                    }
                                }

                                // Parar por tiempo

                                UpdateTimeParal(paralDuraciones[idThread]);
                            }
                        }
                        UpdateTime();
                    }
                }
            }

            // Escribir resultados

            incumbente = *max_element(paralIncumbente.begin(), paralIncumbente.end());
            WriteData(ins, repeticion, (int)accumulate(paralNIter.begin(), paralNIter.end(), 0), nThreads, (int)maxTime);
        }
    }
    return 0;
}
*/
// Single Thread
/*
int main(int argc, char** argv) {

    // Parámetros

    string ins = "";
    double maxTime = 60;
    int iterIni = 1;
    int iterFin = 1;// Inclusivo
    int insIni = 0;
    int insFin = 0; // Inclusivo
    for (int i = 1; i < argc - 1; i += 2) {
        if (argc - 1 >= i + 1) {
            if (string(argv[i]) == "-ins") ins = argv[i + 1];
            else if (string(argv[i]) == "-maxTime") maxTime = atof(argv[i + 1]);
            else if (string(argv[i]) == "-estabilidad") r_estabilidad = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-juntar") r_juntarEspacios = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-ini") iterIni = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-fin") iterFin = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-insIni") insIni = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-insFin") insFin = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-presion") {
                if (string(argv[i + 1]) == "1") r_maxPresionItems = true;
                else r_maxPresionItems = false;
            }
            else if (string(argv[i]) == "-multidrop") r_multidrop = atoi(argv[i + 1]);
            else {
                cout << "Mal en par�metros" << endl;
                return 0;
            }
        }
    }

    // Determinar los clientes

    string ins0 = ins;
    vector<string> misInstancias({ "" });
    if (insFin - insIni >= 1) {
        misInstancias.clear();
        for (int i = insIni; i <= insFin; ++i) {
            misInstancias.push_back(to_string(i));
        }
    }

    // Ciclo de repeticiones

    for (vector<string>::iterator s_it = misInstancias.begin(); s_it < misInstancias.end(); ++s_it) {
        ins = ins0;
        if ((*s_it) != "") {
            ins.replace(ins.find("__"), 2, "_" + *s_it + "_");
        }
        for (int repeticion = iterIni; repeticion <= iterFin; ++repeticion) {
            ReiniciarVariables();

            // Leer datos

            ReadData2(ins);
            int numIter = 0;
            incumbente = Container(C0, false);
            errorArea = 0.5 / (double)(ContenedorDimx * ContenedorDimy);
            double totalPesoAcum = 0;
            int j = 1;
            for (vector<double>::iterator p_it = pesos_acum.begin(); p_it < pesos_acum.end(); ++p_it, ++j) {
                (*p_it) = 1.0 / (double)j;
                totalPesoAcum += *p_it;
            }
            for (vector<double>::iterator p_it = pesos_acum.begin(); p_it < pesos_acum.end(); ++p_it) {
                (*p_it) /= totalPesoAcum;
            }
            for (vector<double>::iterator p_it = pesos_acum.begin() + 1; p_it < pesos_acum.end(); ++p_it) {
                (*p_it) += *(p_it - 1);
            }

            // Total Volumen

            int nCajas = 0;
            nCajas = accumulate(C0.boxest.begin(), C0.boxest.end(), nCajas, [](int const& resp, BoxType const& bt)->int {return resp + bt.cantidad0; });
            totalVolumen = 0;
            for (vector<BoxType>::iterator b_it = C0.boxest.begin(); b_it < C0.boxest.end(); ++b_it) {
                totalVolumen += (*b_it).cantidad0 * (*b_it).volumenReal;
            }
            totalVolumen /= 100.0;
            if (r_multidrop > 0 && totalVolumen > ContenedorVol) { // Se eliminan los clientes que no alcanzan a empacarse por el volumen de los anteriores
                int miVol = 0;
                vector<BoxType>::iterator b_it = C0.boxest.begin();
                int clienteHasta = 0;
                for (; b_it < C0.boxest.end(); ++b_it) {
                    miVol += (*b_it).cantidad0 * (*b_it).volumenReal;
                    if (miVol >= ContenedorVol) {
                        clienteHasta = (*b_it).cliente;
                        break;
                    }
                }
                for (; b_it < C0.boxest.end(); ++b_it) {
                    if ((*b_it).cliente > clienteHasta) {
                        C0.boxest.resize(distance(C0.boxest.begin(), b_it - 1));
                        break;
                    }
                }
            }

            // Adquisición de datos iniciales

            if (r_multidrop > 0) {
                if (r_maxPresionItems) {
                    for (int i = 0; i < alphas.size(); ++i) {
                        ++numIter;
                        Alpha& a = alphas[i];
                        Container C(C0, true);
                        C.alpha = a.val;
                        C.Constructivo_Multidrop_Presion();
                        ActualizarAlpha(a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Multidrop_Presion();
                            C.ConstructivoDeterministico_Multidrop_Presion();
                            ActualizarAlpha(a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Parar por tiempo

                        UpdateTime();
                        if (duracion > maxTime) {
                            break;
                        }
                    }
                }
                else {
                    for (int i = 0; i < alphas.size(); ++i) {
                        ++numIter;
                        Alpha& a = alphas[i];
                        Container C(C0, true);
                        C.alpha = a.val;
                        C.Constructivo_Multidrop();
                        ActualizarAlpha(a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Multidrop();
                            C.ConstructivoDeterministico_Multidrop();
                            ActualizarAlpha(a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Parar por tiempo

                        UpdateTime();
                        if (duracion > maxTime) {
                            break;
                        }
                    }
                }
            }
            else {
                if (r_maxPresionItems) {
                    for (int i = 0; i < alphas.size(); ++i) {
                        ++numIter;
                        Alpha& a = alphas[i];
                        Container C(C0, true);
                        C.alpha = a.val;
                        C.Constructivo_Presion();
                        ActualizarAlpha(a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Presion();
                            C.ConstructivoDeterministico_Presion();
                            ActualizarAlpha(a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Parar por tiempo

                        UpdateTime();
                        if (duracion > maxTime) {
                            break;
                        }
                    }
                }
                else {
                    for (int i = 0; i < alphas.size(); ++i) {
                        ++numIter;
                        Alpha& a = alphas[i];
                        Container C(C0, true);
                        C.alpha = a.val;
                        C.Constructivo();
                        ActualizarAlpha(a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas();
                            C.ConstructivoDeterministico();
                            ActualizarAlpha(a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Parar por tiempo

                        UpdateTime();
                        if (duracion > maxTime) {
                            break;
                        }
                    }
                }
            }

            // Resto de la metaheurística

            if (duracion > maxTime || !incumbente.hayCajasPorEmpacar) {
                WriteData(ins, repeticion, numIter, 1, 0);
                continue;
            }
            if (r_multidrop > 0) {
                if (r_maxPresionItems) {
                    while (duracion < maxTime) {
                        ++numIter;
                        vector<Alpha>::iterator a = DeterminarAlpha();
                        Container C = Container(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo_Multidrop_Presion();
                        ActualizarAlpha(*a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Multidrop_Presion();
                            C.ConstructivoDeterministico_Multidrop_Presion();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Actualizar tiempo

                        UpdateTime();
                    }
                }
                else {
                    while (duracion < maxTime) {
                        ++numIter;
                        vector<Alpha>::iterator a = DeterminarAlpha();
                        Container C = Container(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo_Multidrop();
                        ActualizarAlpha(*a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Multidrop();
                            C.ConstructivoDeterministico_Multidrop();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Actualizar tiempo

                        UpdateTime();
                    }
                }
            }
            else {
                if (r_maxPresionItems) {
                    while (duracion < maxTime) {
                        ++numIter;
                        vector<Alpha>::iterator a = DeterminarAlpha();
                        Container C = Container(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo_Presion();
                        ActualizarAlpha(*a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Presion();
                            C.ConstructivoDeterministico_Presion();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Parar por tiempo

                        UpdateTime();
                    }
                }
                else {
                    while (duracion < maxTime) {
                        ++numIter;
                        vector<Alpha>::iterator a = DeterminarAlpha();
                        Container C = Container(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo();
                        ActualizarAlpha(*a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas();
                            C.ConstructivoDeterministico();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Actualizar tiempo

                        UpdateTime();
                    }
                }
            }

            // Escribir resultados

            //cout << "total iteraciones = " << to_string(holi) << endl;
            WriteData(ins, repeticion, numIter, 1, 0);
        }
    }
    return 0;
}
*/
// Single Thread MLBR
/*
int main(int argc, char** argv) {

    // Parámetros

    string ins = "";
    int iterIni = 1;
    int iterFin = 1;// Inclusivo
    int insIni = 0;
    int insFin = 0; // Inclusivo
    int maxIter = 50000;
    for (int i = 1; i < argc - 1; i += 2) {
        if (argc - 1 >= i + 1) {
            if (string(argv[i]) == "-ins") ins = argv[i + 1];
            else if (string(argv[i]) == "-maxIter") maxIter = atof(argv[i + 1]);
            else if (string(argv[i]) == "-estabilidad") r_estabilidad = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-juntar") r_juntarEspacios = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-ini") iterIni = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-fin") iterFin = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-insIni") insIni = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-insFin") insFin = atoi(argv[i + 1]);
            else if (string(argv[i]) == "-presion") {
                if (string(argv[i + 1]) == "1") r_maxPresionItems = true;
                else r_maxPresionItems = false;
            }
            else if (string(argv[i]) == "-multidrop") r_multidrop = atoi(argv[i + 1]);
            else {
                cout << "Mal en par�metros" << endl;
                return 0;
            }
        }
    }

    // Determinar los clientes

    string ins0 = ins;
    vector<string> misInstancias({ "" });
    if (insFin - insIni >= 1) {
        misInstancias.clear();
        for (int i = insIni; i <= insFin; ++i) {
            misInstancias.push_back(to_string(i));
        }
    }

    // Ciclo de repeticiones

    for (vector<string>::iterator s_it = misInstancias.begin(); s_it < misInstancias.end(); ++s_it) {
        ins = ins0;
        if ((*s_it) != "") {
            ins.replace(ins.find("__"), 2, "_" + *s_it + "_");
        }
        for (int repeticion = iterIni; repeticion <= iterFin; ++repeticion) {
            ReiniciarVariables();

            // Leer datos

            ReadData2(ins);
            int numIter = 0;
            incumbente = Container(C0, false);
            errorArea = 0.5 / (double)(ContenedorDimx * ContenedorDimy);
            double totalPesoAcum = 0;
            int j = 1;
            for (vector<double>::iterator p_it = pesos_acum.begin(); p_it < pesos_acum.end(); ++p_it, ++j) {
                (*p_it) = 1.0 / (double)j;
                totalPesoAcum += *p_it;
            }
            for (vector<double>::iterator p_it = pesos_acum.begin(); p_it < pesos_acum.end(); ++p_it) {
                (*p_it) /= totalPesoAcum;
            }
            for (vector<double>::iterator p_it = pesos_acum.begin() + 1; p_it < pesos_acum.end(); ++p_it) {
                (*p_it) += *(p_it - 1);
            }

            // Total Volumen

            int nCajas = 0;
            nCajas = accumulate(C0.boxest.begin(), C0.boxest.end(), nCajas, [](int const& resp, BoxType const& bt)->int {return resp + bt.cantidad0; });
            totalVolumen = 0;
            for (vector<BoxType>::iterator b_it = C0.boxest.begin(); b_it < C0.boxest.end(); ++b_it) {
                totalVolumen += (*b_it).cantidad0 * (*b_it).volumenReal;
            }
            totalVolumen /= 100.0;
            if (r_multidrop > 0 && totalVolumen > ContenedorVol) { // Se eliminan los clientes que no alcanzan a empacarse por el volumen de los anteriores
                int miVol = 0;
                vector<BoxType>::iterator b_it = C0.boxest.begin();
                int clienteHasta = 0;
                for (; b_it < C0.boxest.end(); ++b_it) {
                    miVol += (*b_it).cantidad0 * (*b_it).volumenReal;
                    if (miVol >= ContenedorVol) {
                        clienteHasta = (*b_it).cliente;
                        break;
                    }
                }
                for (; b_it < C0.boxest.end(); ++b_it) {
                    if ((*b_it).cliente > clienteHasta) {
                        C0.boxest.resize(distance(C0.boxest.begin(), b_it - 1));
                        break;
                    }
                }
            }

            // Adquisición de datos iniciales

            if (r_multidrop > 0) {
                if (r_maxPresionItems) {
                    for (int i = 0; i < alphas.size(); ++i) {
                        ++numIter;
                        Alpha& a = alphas[i];
                        Container C(C0, true);
                        C.alpha = a.val;
                        C.Constructivo_Multidrop_Presion();
                        ActualizarAlpha(a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Multidrop_Presion();
                            C.ConstructivoDeterministico_Multidrop_Presion();
                            ActualizarAlpha(a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Parar por iter

                        if (numIter >= maxIter) {
                            break;
                        }
                    }
                }
                else {
                    for (int i = 0; i < alphas.size(); ++i) {
                        ++numIter;
                        Alpha& a = alphas[i];
                        Container C(C0, true);
                        C.alpha = a.val;
                        C.Constructivo_Multidrop();
                        ActualizarAlpha(a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Multidrop();
                            C.ConstructivoDeterministico_Multidrop();
                            ActualizarAlpha(a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Parar por iter

                        if (numIter >= maxIter) {
                            break;
                        }
                    }
                }
            }
            else {
                if (r_maxPresionItems) {
                    for (int i = 0; i < alphas.size(); ++i) {
                        ++numIter;
                        Alpha& a = alphas[i];
                        Container C(C0, true);
                        C.alpha = a.val;
                        C.Constructivo_Presion();
                        ActualizarAlpha(a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Presion();
                            C.ConstructivoDeterministico_Presion();
                            ActualizarAlpha(a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Parar por iter

                        if (numIter >= maxIter) {
                            break;
                        }
                    }
                }
                else {
                    for (int i = 0; i < alphas.size(); ++i) {
                        ++numIter;
                        Alpha& a = alphas[i];
                        Container C(C0, true);
                        C.alpha = a.val;
                        C.Constructivo();
                        ActualizarAlpha(a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas();
                            C.ConstructivoDeterministico();
                            ActualizarAlpha(a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }

                        // Parar por iter

                        if (numIter >= maxIter) {
                            break;
                        }
                    }
                }
            }

            // Resto de la metaheurística

            if (numIter >= maxIter || !incumbente.hayCajasPorEmpacar) {
                WriteData(ins, repeticion, numIter, 1, 0);
                continue;
            }
            if (r_multidrop > 0) {
                if (r_maxPresionItems) {
                    while (numIter < maxIter) {
                        ++numIter;
                        vector<Alpha>::iterator a = DeterminarAlpha();
                        Container C = Container(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo_Multidrop_Presion();
                        ActualizarAlpha(*a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Multidrop_Presion();
                            C.ConstructivoDeterministico_Multidrop_Presion();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }
                    }
                }
                else {
                    while (numIter < maxIter) {
                        ++numIter;
                        vector<Alpha>::iterator a = DeterminarAlpha();
                        Container C = Container(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo_Multidrop();
                        ActualizarAlpha(*a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Multidrop();
                            C.ConstructivoDeterministico_Multidrop();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }
                    }
                }
            }
            else {
                if (r_maxPresionItems) {
                    while (numIter < maxIter) {
                        ++numIter;
                        vector<Alpha>::iterator a = DeterminarAlpha();
                        Container C = Container(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo_Presion();
                        ActualizarAlpha(*a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas_Presion();
                            C.ConstructivoDeterministico_Presion();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }
                    }
                }
                else {
                    while (numIter < maxIter) {
                        ++numIter;
                        vector<Alpha>::iterator a = DeterminarAlpha();
                        Container C = Container(C0, true);
                        C.alpha = (*a).val;
                        C.Constructivo();
                        ActualizarAlpha(*a, C.utilizacion);
                        int valorActualizacion = ActualizarIncumbente(C);
                        if (valorActualizacion == -1) {
                            break;
                        }
                        if (valorActualizacion == 1 || C.utilizacion > threshold_val) {

                            // Búsqueda local

                            C.RemoverPiezas();
                            C.ConstructivoDeterministico();
                            ActualizarAlpha(*a, C.utilizacion);
                            if (ActualizarIncumbente(C) == -1) {
                                break;
                            }
                            threshold_val = incumbente.utilizacion;
                            threshold_Iter = 0;
                        }
                        else {
                            ++threshold_Iter;
                            if (threshold_Iter == threshold_MaxIter) {
                                threshold_val *= 0.8;
                                threshold_Iter = 0;
                            }
                        }
                    }
                }
            }

            // Escribir resultados

            //cout << "total iteraciones = " << to_string(holi) << endl;
            WriteData(ins, repeticion, numIter, 1, 0);
        }
    }
    return 0;
}
*/
//END