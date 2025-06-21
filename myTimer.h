#ifndef MYTIMER_H
#define MYTIMER_H

#include <chrono>

/*!
 * \brief Timer - всё необходимое для измерения времени на ЦП
 */
namespace Timer {

/// Далее объявлено всё, что касается измерения времени
/*!
 * \brief time_var - тип данных, хранящий временную метку
 */
typedef std::chrono::high_resolution_clock::time_point time_var;

#define time_now() ( std::chrono::high_resolution_clock::now() ); /// текущий момент времени
#define time_taken(p1, p2) ( std::chrono::duration_cast<std::chrono::nanoseconds>( (p2) - (p1) ).count() / 1e9 ); /// подсчёт длительности

#define prec std::setprecision(std::numeric_limits<mpComplex_t::value_type>::digits10)

};

#endif // MYTIMER_H
