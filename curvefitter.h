#ifndef CURVEFITTER_H
#define CURVEFITTER_H

#include <cmath>
#include <cstring>
#include <vector>
#include  "mpfit.h"

typedef double (*custom_calc_func)(const double x,double coef[]);
static int mpfit_routine(int m, int n, double *p, double *dy, double **dvec, void *vars); // mp_func  used in mpfit() as callback;


struct vars_struct {
    double *x;
    double *y;
    double *ey;   
    custom_calc_func func;
};


class CurveFitter
{
public:

    enum class FittingStatus
    {
        FS_OK=0,
        FS_SOURCE_DATA_NOT_SET,
        FS_SOURCE_NOT_ENOUGH_ELEMENTS,
        FS_SOURCE_CANT_BE_NEGATIVE,
        FS_SOURCE_VAL_CANT_BE_ZERO,
        FS_SOURCE_VAL_MUST_BE_POSITIVE,
        FS_CUSTOM_FUNC_NOT_SET,
        FS_COEF_NUMB_TO_BIG,
        FS_RESULT_Y_NULL,
        FS_INPUT_X_NULL,
        FS_SOURCE_X_DATA_NOT_ORDERED
    };

    enum class FittingModel
    {
        LINEAR=0,
        POLYNOMIAL2,
        POLYNOMIAL3,
        POLYNOMIAL4,
        POLYNOMIAL5,
        LOGREGRESSION,
        EXPREGRESSION,
        PIECWEISELINEAR,
        USERDEFINED
    };

    CurveFitter();
    CurveFitter(const double x_source[],const double y_source[], size_t element_number);
    CurveFitter(const CurveFitter&) = delete;
    CurveFitter(const CurveFitter&&) = delete;
    ~CurveFitter()=default;

    CurveFitter& operator= (const  CurveFitter&) =delete;
    CurveFitter& operator= (const  CurveFitter&&) =delete;

    FittingStatus set_source_data_points(const double x_source[],const double y_source[], size_t element_number);
    FittingStatus calculate(const double x[],size_t points_to_calc, FittingModel model,double result_y[]);
    FittingStatus set_custom_fitting_model(custom_calc_func func, int coef_numb);
    FittingStatus calculate_coef(FittingModel model, std::vector<double>& coef );


    void realese_source_data(void);
    void release_cutom_funct(void);

    const double* get_source_x(void);
    const double* get_source_y(void);
    size_t get_xy_elem_numb(void);



private:

    FittingStatus validate_data(const double x[],double result_y[]);
    FittingStatus validate_data (FittingModel model);
    FittingStatus validate_source_data (void);

    double calc_piecweise_point(double x_point);
    void SetRoutineParam(FittingModel mode);

const double *x_ptr, *y_ptr;
size_t element_number;

mp_config config;

vars_struct param;

mp_result result;

int custom_func_coef_numb, routine_coef_numb;

custom_calc_func routine,custom_func;

const static int MAX_PARAM_NUMB=20;
double param_arr[MAX_PARAM_NUMB];


};

#endif // CURVEFITTER_H
