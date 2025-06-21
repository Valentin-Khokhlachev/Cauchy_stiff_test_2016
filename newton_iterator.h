#ifndef NEWTON_ITERATOR_H
#define NEWTON_ITERATOR_H

#include <QCoreApplication>
#include <QtCore>
#include <QVector>
#include <QDebug>
#include <QString>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <functional>

#include <eigen3/Eigen/Eigen>

#include "mpTypes.h"
#include "myTimer.h"
#include "csvstream.h"
#include "config_parser.h"


using namespace std;
using namespace Eigen;
using namespace mpTypes;

namespace Newton {

class Newton_iterator
{
public:
    Newton_iterator();

    void setPrecForPrint(int prc);
    int getPrecForPrint();

    void setEps(string realPart, string imagPart);
    mpComplex_t getEps();

    void setTheta(string realPart, string imagPart);
    mpComplex_t getTheta();

    void setMaxIterations(int maxIt);
    int getMaxIterations();

    void setAlpha(string realPart, string imagPart);
    mpComplex_t getAlpha();

    /*!
     * \brief Newton_method - метод, реализующий Ньютоновский итерационный процесс
     * \param f - функция, корень которой ищем ( \vec{f}(\vec{x})= 0 )
     * \param J_f - матрица Якоби для \vec{f}(\vec{x}):
     * (вектор-столбец, составленный из градиентов-строк)
     * \param v_prev - начальное приближение \vec{x}_{0}
     * \return приближенное решение \vec{\tilde{x}}
     */
    VectorX<mpComplex_t> Newton_method(function<VectorX<mpComplex_t>(VectorX<mpComplex_t>)> f,
                                        function<MatrixX<mpComplex_t>(VectorX<mpComplex_t>)> J_f,
                                        const VectorX<mpComplex_t> &v_prev);

private:
    /*!
     * \brief precForPrint - точность, с которой делается отладочный вывод в данном модуле
     */
    unsigned int precForPrint = 8u;

    /*!
     * \brief eps - параметр для условия выхода из итерационного Ньютоновского процесса
     */
    mpComplex_t eps {"1e-16", "0.0"}; /// = 1e-16

    /*!
     * \brief theta - параметр из обобщённого метода Ньютона
     */
    mpComplex_t theta {"1.0", "0.0"};

    /*!
     * \brief maxIterations - максимальное количество итераций по обобщённому методу Ньютона
     */
    int maxIterations = 120;

    /*!
     * \brief alpha - множитель для случайной части
     * начального приближения для метода Ньютона
     * (если 0, то начальное приближение строго детерминированное)
     */
    mpComplex_t alpha {"0.0", "0.0"};
};

};

#endif // NEWTON_ITERATOR_H
