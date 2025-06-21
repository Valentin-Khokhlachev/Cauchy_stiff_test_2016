#ifndef IVP_CAUCHY_SOLVER_H
#define IVP_CAUCHY_SOLVER_H

#include <QCoreApplication>
#include <QtCore>
#include <QVector>
#include <QDebug>
#include <QString>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <functional>

#include <eigen3/Eigen/Eigen>

#include "mpTypes.h"
#include "myTimer.h"
#include "newton_iterator.h"
#include "csvstream.h"
#include "config_parser.h"



/*!
 * \brief Cauchy - всё, что касается решения задачи Коши для скалярного уравнения
 * в аргументе длина дуги при помощи экспоненциальной разностной схемы
 */
namespace Cauchy {

using namespace std;
using namespace Eigen;
using namespace mpTypes;
// using boost::multiprecision::cpp_dec_float;
// using boost::multiprecision::cpp_complex;

class ivpCauchySolver
{
public:
    /*!
     * \brief ivp_Cauchy_solver - конструктор
     */
    ivpCauchySolver();

    /*!
     * \brief ~ivp_Cauchy_solver - деструктор
     */
    ~ivpCauchySolver();

    /*!
     * \brief solve - метод, рещающий вшитую в класс начальную задачу
     */
    void solve();

    /*!
     * \brief write_full_dump - вывод всех данных
     * в csv файл "vpa_stiff_test_2016_arg_arc.csv"
     */
    void writeData(const int pr);

private:
    /*!
     * \brief PI - число пи с раширенной разрядностью
     */
    const mpFloat_t PI = mpFloat_t("4.0") * atan( mpFloat_t("1.0") );

    /*!
     * \brief PI_С - число пи с раширенной разрядностью
     */
    const mpComplex_t PI_C {PI, "0.0"};

    /*!
     * \brief I - мнимая единица
     */
    const mpComplex_t I {"0.0", "1.0"};

    /*!
     * \brief N - число точек в каждой из сеток
     * (число точек в сетке на вещественной прямой совпадает с
     *  числом точек по контуру в комплексной плоскости)
     */
    int N = 8;

    /*!
     * \brief l_0 - начальная точка по длине дуги
     */
    mpComplex_t l_0 {"0.0", "0.0"};

    /*!
     * \brief L - конечная точка по длине дуги
     */
    mpComplex_t L {"3.0", "0.0"};

    /*!
     * \brief t_0 - начальная точка по времени
     */
    mpComplex_t t_0 {"0.0", "0.0"};

    /*!
     * \brief alpha - параметр для генерации начального приближения в методе Ньтона
     * (если равно нулю, то первое приближение полностью детерминированная величина)
     */
    mpComplex_t alpha {"0.0", "0.0"};

    /*!
     * \brief uInit - вектор из начальных значений
     */
    VectorX<mpComplex_t> uInit;

    /*!
     * \brief pA - параметр a из исходной задачи
     */
    mpComplex_t pA {"1.0", "0.0"};

    /*!
     * \brief pLambda0 - параметр lambda_0 из исходной задачи
     */
    mpComplex_t pLambda0 {"10.0", "0.0"};

    /*!
     * \brief rho - радиус окружности ( должен быть > 0.5*(L-l_0) <- для длины дуги  )
     */
    mpComplex_t rho = mpFloat_t("1.5") * (L - l_0);

    /*!
     * \brief l - равномерная сетка по длине дуги интегральной кривой
     */
    VectorX<mpComplex_t> l;

    /*!
     * \brief h_phi - шаг равномерной сетки по углу в полярной системе координат
     */
    mpComplex_t h_phi = mpFloat_t("2.0") * PI_C / mpFloat_t(N - 1);

    /*!
     * \brief phi - равномерная сетка по углу в полярной системе координат
     */
    VectorX<mpComplex_t> phi;

    /*!
     * \brief z - сетка на замкнутом контуре в комплексной плоскости
     */
    VectorX<mpComplex_t> z;

    /*!
     * \brief h_z - шаги, приписанне к узлам сетки z
     */
    VectorX<mpComplex_t> h_z;

    /*!
     * \brief uWaveInit - начальное приближение сеточного аналитического
     * продолжения сеточной функции u для метода Ньютона
     */
    VectorX<mpComplex_t> uWaveInit;

    /*!
     * \brief uWave - сеточное аналитическое продолжение сеточной функции u
     */
    VectorX<mpComplex_t> uWave;

    /*!
     * \brief u - искомая сеточная функция
     */
    VectorX<mpComplex_t> u;

    /*!
     * \brief uExact - точный ответ
     */
    VectorX<mpComplex_t> uExact;




    /*!
     * \brief file_path - относительный путь к папке, где будет лежать файл
     */
    string filePath = "../../results/";

    /*!
     * \brief file_name - имя файл, куда будем выводить результаты
     */
    string fileName = "vpa_stiff_test_2016_arg_arc.csv";

    /*!
     * \brief file - относительный путь к файлу
     */
    string file = filePath + fileName;

    /*!
     * \brief ofs - поток вывода в файл (из стандартной библиотеки)
     */
    std::ofstream ofs;

    /*!
     * \brief csvfile - экземпляр класса-csv-файла
     */
    // unique_ptr<csvstream> csvfile;



    /*!
     * \brief iter - можкль, осуществляющий Ньютоновский итерационный процесс
     */
    unique_ptr<Newton::Newton_iterator> iter;



    /*!
     * \brief contour - генерация сетки в комплексной плоскости
     * \param t_s - начало отрезка на вещественной оси
     * \param t_f - конец отрезка на вещественной оси
     * \param rho - радиус окружности
     * \param phi - сетка по углу (вектор)
     * \return сетку в комплексной плоскости (вектор)
     */
    VectorX<mpComplex_t> contour(const mpComplex_t t_s,
                                 const mpComplex_t t_f,
                                 const mpComplex_t rho,
                                 const Ref<VectorX<mpComplex_t>> phi);

    /*!
     * \brief h_contour - генерация шагов, приписанных к узлам сетки (вектор)
     * \param rho - радиус окружности
     * \param phi - равномерная сетка по углу (вектор)
     * \param h_phi - шаг сетки по углу
     * \return вектор из шагов, приписанных к узлам сетки в комплексной плоскости
     */
    VectorX<mpComplex_t> hContour(const mpComplex_t rho,
                                   const VectorX<mpComplex_t> phi,
                                   const mpComplex_t h_phi);

    /*!
     * \brief mstr - переводит mpComplex_t в строку, понятную matlab'у
     * \param z - экземпляр типа mpComplex_t
     * \return строка, отформатированная для чтения matlab'ом
     */
    string mstr(const mpComplex_t z, const int pr);

    /*!
     * \brief Newton_init - функция, генерирующая начальное приближение
     * для метода Ньютона
     * \param x_s - стартовая точка на вещественной прямой
     * \param x_f - финишная точка на вещественной прямой
     * \param N - кол-во точек (у нас = кол-ву шагов (замкнутость контура))
     * \return начальное приближение для метода Ньютона
     */
    void NewtonInit();

    /*!
     * \brief recalc - метод восстанавливающий искомую сеточную функцию
     * по её аналитическому продолжению
     * \return
     */
    void reCalc();

    /// Правая часть - первая строка
    mpComplex_t rightHand0(VectorX<mpComplex_t> x);
    /// Правая часть - второая строка
    mpComplex_t rightHand1(VectorX<mpComplex_t> x);
    /// Частная производная первой строки правой части по первому аргументу
    mpComplex_t rightHandDerivative00(VectorX<mpComplex_t> x);
    /// Частная производная первой строки правой части по второму аргументу
    mpComplex_t rightHandDerivative01(VectorX<mpComplex_t> x);
    /// Частная производная второй строки правой части по первому аргументу
    mpComplex_t rightHandDerivative10(VectorX<mpComplex_t> x);
    /// Частная производная второй строки правой части по второму аргументу
    mpComplex_t rightHandDerivative11(VectorX<mpComplex_t> x);

    /*!
     * \brief F - вектор-функция правая часть
     * \param x - аргумент
     * \return значение вектор-функции F
     */
    VectorX<mpComplex_t> F (VectorX<mpComplex_t> x);

    /*!
     * \brief JF - матрица Якоби для вектор-функции F
     * \param x - аргумент
     * \return матрица Якоби для вектор-функции F
     */
    MatrixX<mpComplex_t> JF (VectorX<mpComplex_t> x);

    /*!
     * \brief exactSolution - точное решение начальной задачи
     * \param x - аргумент
     */
    void exactSolution (VectorX<mpComplex_t> x);

    /*!
     * \brief relErr - метод, вычисляющий относительную погрешность
     * приближенного решения в норме L2
     * \param x - вектор из аргументов - времен
     * \return относительная погрешность в норме L2
     */
    mpComplex_t relErr ();

    /*!
     * \brief testFunc - некоторая вектор-функция (для теста)
     * \param x - аргумент
     * \return значение вектор-функции testFunc
     */
    VectorX<mpComplex_t> testFunc (VectorX<mpComplex_t> x);

    /*!
     * \brief testJFunc - матрица Якоби для вектор-функции testFunc
     * \param x - аргумент
     * \return матрица Якоби для вектор-функции testFunc
     */
    MatrixX<mpComplex_t> testJFunc (VectorX<mpComplex_t> x);
};

}

#endif // IVP_CAUCHY_SOLVER_H
