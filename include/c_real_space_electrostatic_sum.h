#ifndef __C_REAL_SPACE_ELECTROSTATIC_SUM_H__
#define __C_REAL_SPACE_ELECTROSTATIC_SUM_H__

extern "C"
void c_real_space_electrostatic_sum_energy(
        const double* a1, const double* a2, const double* a3,
        const int* num,
        const double* rx, const double* ry, const double* rz,
        const double* z,
        const double* rc,
        const double* rd,
        double* e);

extern "C"
void c_real_space_electrostatic_sum_force( 
        const double* a1, const double* a2, const double* a3,
        const int* num,
        const double* rx, const double* ry, const double* rz,
        const double* z,
        const double* rc,
        const double* rd,
        double* fx, double* fy, double* fz);

extern "C"
void c_real_space_electrostatic_sum_stress(
        const double* a1, const double* a2, const double* a3,
        const int* num,
        const double* rx, const double* ry, const double* rz,
        const double* z,
        const double* rc,
        const double* rd,
        double* s);

#endif // __C_REAL_SPACE_ELECTROSTATIC_SUM_H
