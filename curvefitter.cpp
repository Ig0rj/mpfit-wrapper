#include "curvefitter.h"

CurveFitter::CurveFitter()
{
    this->x_ptr=nullptr;
    this->y_ptr=nullptr;
    this->element_number=0;
    param.x=nullptr;
    param.y=nullptr;
    param.ey=nullptr;
    param.func=nullptr;

    this->custom_func_coef_numb=0;
    this->routine_coef_numb=0;
}

CurveFitter::CurveFitter(const double x_source[], const double y_source[], size_t element_number) :CurveFitter()
{
    this->x_ptr=x_source;
    this->y_ptr=y_source;
    this->element_number=element_number;

    this->param.x=const_cast<double*>(x_source);
    this->param.y=const_cast<double*>(y_source);
}


 CurveFitter::~CurveFitter()
{



}

const double* CurveFitter::get_source_x(void)
{
        return this->x_ptr;
}

const double* CurveFitter::get_source_y(void)
{
    return this->y_ptr;
}

size_t CurveFitter::get_xy_elem_numb(void)
{
    return this->element_number;
}

void  CurveFitter::realese_source_data(void)
{
    this->x_ptr=nullptr;
    this->y_ptr=nullptr;
    this->element_number=0;
    param.x=nullptr;
    param.y=nullptr;
    param.ey=nullptr;

}
void  CurveFitter::release_cutom_funct(void)
{
    this->custom_func_coef_numb=0;
    this->routine_coef_numb=0;
}

 CurveFitter::FittingStatus CurveFitter::validate_data(const double x[],double result_y[])
{

     if(result_y==nullptr)
         return FittingStatus::FS_RESULT_Y_NULL;

     if(x==nullptr)
         return   FittingStatus::FS_INPUT_X_NULL;


    return FittingStatus::FS_OK;
}

CurveFitter::FittingStatus CurveFitter::validate_source_data(void)
{
    if ((this->x_ptr==nullptr) || (this->y_ptr==nullptr) || (this->element_number==0))
        return  FittingStatus::FS_SOURCE_DATA_NOT_SET;

    if(x_ptr[0]<x_ptr[1])
    {
        for(auto i=0;i<this->element_number-1;i++)
        {
            if(x_ptr[i]>=x_ptr[i+1])
                return FittingStatus::FS_SOURCE_X_DATA_NOT_ORDERED;
        }
    }
    else
    {
        for(auto i=0;i<this->element_number-1;i++)
        {
            if(x_ptr[i]<=x_ptr[i+1])
                return FittingStatus::FS_SOURCE_X_DATA_NOT_ORDERED;
        }
    }

     return FittingStatus::FS_OK;
}

CurveFitter::FittingStatus  CurveFitter::validate_data(FittingModel model)
{
    static double zero_val;

    CurveFitter::FittingStatus res=validate_source_data();
    if(res!= FittingStatus::FS_OK)
        return res;

    switch(model)
    {

    case FittingModel::LINEAR:

        if(this->element_number<2)
            return FittingStatus::FS_SOURCE_NOT_ENOUGH_ELEMENTS;

        break;

    case FittingModel::POLYNOMIAL2:

        if(this->element_number<3)
            return FittingStatus::FS_SOURCE_NOT_ENOUGH_ELEMENTS;

        break;

    case FittingModel::POLYNOMIAL3:

        if(this->element_number<4)
            return FittingStatus::FS_SOURCE_NOT_ENOUGH_ELEMENTS;

        break;

    case FittingModel::POLYNOMIAL4:

        if(this->element_number<5)
            return FittingStatus::FS_SOURCE_NOT_ENOUGH_ELEMENTS;

        break;

    case FittingModel::POLYNOMIAL5:

        if(this->element_number<6)
            return FittingStatus::FS_SOURCE_NOT_ENOUGH_ELEMENTS;

        break;

    case FittingModel::LOGREGRESSION:

        if(this->element_number<2)
            return FittingStatus::FS_SOURCE_NOT_ENOUGH_ELEMENTS;

        for(int i=0;i<this->element_number;i++)
        {
            if(this->x_ptr[i]<=zero_val)
                return FittingStatus::FS_SOURCE_VAL_MUST_BE_POSITIVE;
        }

        break;

    case FittingModel::EXPREGRESSION:

        if(this->element_number<2)
            return FittingStatus::FS_SOURCE_NOT_ENOUGH_ELEMENTS;

        for(int i=0;i<this->element_number;i++)
        {
            if(this->x_ptr[i]<=zero_val)
                return FittingStatus::FS_SOURCE_VAL_MUST_BE_POSITIVE;
        }

        break;


    case FittingModel::PIECWEISELINEAR:

        if(this->element_number<2)
            return FittingStatus::FS_SOURCE_NOT_ENOUGH_ELEMENTS;

        break;

    case FittingModel::USERDEFINED:

        if(this->custom_func==nullptr)
            return FittingStatus:: FS_CUSTOM_FUNC_NOT_SET;

        if(this->custom_func_coef_numb>MAX_PARAM_NUMB)
            return FittingStatus:: FS_COEF_NUMB_TO_BIG;

        break;

    }

    return FittingStatus::FS_OK;
}



CurveFitter::FittingStatus  CurveFitter::set_custom_fitting_model(custom_calc_func func, int coef_numb)
{

    if(coef_numb>MAX_PARAM_NUMB)
       return FittingStatus:: FS_COEF_NUMB_TO_BIG;

    this->custom_func=func;
    this->custom_func_coef_numb=coef_numb;

    return FittingStatus::FS_OK;

}
CurveFitter::FittingStatus CurveFitter::set_source_data_points(const double x_source[],const double y_source[], size_t element_number)
{

    if(element_number<2)
        return FittingStatus::FS_SOURCE_NOT_ENOUGH_ELEMENTS;

    this->x_ptr=const_cast<double*>(x_source);
    this->y_ptr=const_cast<double*>(y_source);
    this->element_number=element_number;

    this->param.x=const_cast<double*>(x_source);
    this->param.y=const_cast<double*>(y_source);

    return FittingStatus::FS_OK;
}

void CurveFitter::SetRoutineParam(FittingModel model)
{

    std::memset(&result,0,sizeof(mp_result));

    std::memset(&config, 0, sizeof(config));
    config.maxiter = 1000;


    for(int i=0;i<MAX_PARAM_NUMB;i++)
        param_arr[i]=1.0;


    switch(model)
    {

    case FittingModel::LINEAR:

        this->routine =[](const double x,double p[])->double{ return p[0] + p[1]*x;};
        this-> routine_coef_numb=2;
        break;

    case FittingModel::POLYNOMIAL2:

        this->routine=[](const double x,double p[])->double{ return (p[0]*pow(x,2))+(p[1]*x)+p[2];};
        this->routine_coef_numb=3;
        break;

    case FittingModel::POLYNOMIAL3:

        this->routine=routine=[](const double x,double p[])->double{ return (p[0]*pow(x,3))+(p[1]*pow(x,2))+(p[2]*x)+p[3];};
        this->routine_coef_numb=4;
        break;

    case FittingModel::POLYNOMIAL4:

        this->routine=[](const double x,double p[])->double{return (p[0]*pow(x,4))+(p[1]*pow(x,3))+(p[2]*pow(x,2))+(p[3]*x)+p[4];};
        this->routine_coef_numb=5;

        break;

    case FittingModel::POLYNOMIAL5:

        this->routine=[](const double x,double p[])->double{ return (p[0]*pow(x,5))+(p[1]*pow(x,4))+(p[2]*pow(x,3))+(p[3]*pow(x,2))+(p[4]*x)+p[5];};
        this->routine_coef_numb=6;
        break;

    case FittingModel::LOGREGRESSION:

        this->routine =[](const double x,double p[])->double{return p[0]*log(x)+p[1];};
        this->routine_coef_numb=2;
        break;

    case FittingModel::EXPREGRESSION:

        this->routine =[](const double x,double p[])->double{return p[0]*pow(x,p[1]);};
        this->routine_coef_numb=2;
        break;


    case FittingModel::PIECWEISELINEAR:

        this->routine =[](const double x,double p[])->double{ return p[0] + p[1]*x;};
        this-> routine_coef_numb=2;
        break;

    case FittingModel::USERDEFINED:

        this->routine=custom_func;
        this->routine_coef_numb=this->custom_func_coef_numb;
        break;

    }
}


CurveFitter::FittingStatus CurveFitter::calculate_coef(FittingModel model, std::vector<double>& coef )
{

    FittingStatus ret=validate_data(model);
    if(ret!= FittingStatus::FS_OK)
        return ret;

    SetRoutineParam(model);

    double* ey= new double[this->element_number];

    for (int i=0; i<this->element_number; i++)
        ey[i] = 0.07;


    param.ey = ey;
    param.func=this->routine;

    if(model!=FittingModel::PIECWEISELINEAR)
    {
            mpfit(mpfit_routine,this->element_number,  this->routine_coef_numb, param_arr, 0, &config, (void *) &param, &result);
            coef.clear();
            coef.reserve(this-> routine_coef_numb);
            for(auto i=0;i<this->routine_coef_numb;i++)
                coef.push_back(param_arr[i]);
    }
    else
    {
        coef.clear();
    }

    delete[] ey;

    return FittingStatus::FS_OK;
}


CurveFitter::FittingStatus CurveFitter::calculate(const double x[],size_t points_to_calc, FittingModel model,double result_y[])
{
     FittingStatus ret=validate_data(x,result_y);

     if(ret!= FittingStatus::FS_OK)
         return ret;

     ret=validate_data(model);
     if(ret!= FittingStatus::FS_OK)
        return ret;

   SetRoutineParam(model);

    double* ey= new double[this->element_number];

    for (int i=0; i<this->element_number; i++)
        ey[i] = 0.07;


    param.ey = ey;
    param.func=this->routine;

    if(model!=FittingModel::PIECWEISELINEAR)
    {

     mpfit(mpfit_routine,this->element_number,  this->routine_coef_numb, param_arr, 0, &config, (void *) &param, &result);
        for(auto i=0;i<points_to_calc;i++)
            result_y[i]=routine(x[i],param_arr);
    }
    else
    {
        for(auto i=0;i<points_to_calc;i++)
            result_y[i]=calc_piecweise_point(x[i]);
    }

    delete[] ey;

    return FittingStatus::FS_OK;

}

double CurveFitter::calc_piecweise_point( double x_point)
{

    if(x_ptr[0]<x_ptr[1])
    {

        if(x_point<=this->x_ptr[1])
        {
            param.x=(double*)this->x_ptr;
            param.y=(double*)this->y_ptr;
        }

        if(x_point>=this->x_ptr[this->element_number-2])
        {
            param.x=(double*)&this->x_ptr[this->element_number-2];
            param.y=(double*)&this->y_ptr[this->element_number-2];
        }

        if((x_point>this->x_ptr[1]) && (x_point<this->x_ptr[this->element_number-2]))
        {
            for(int i=1;i<this->element_number-2;i++)
            {
                if((x_point>=this->x_ptr[i]) && (x_point<=this->x_ptr[i+1]))
                {
                    param.x=(double*)&this->x_ptr[i];
                    param.y=(double*)&this->y_ptr[i];
                }
            }
        }
    }
    else
    {
        if(x_point>=this->x_ptr[1])
        {
            param.x=(double*)this->x_ptr;
            param.y=(double*)this->y_ptr;
        }

        if(x_point<=this->x_ptr[this->element_number-2])
        {
            param.x=(double*)&this->x_ptr[this->element_number-2];
            param.y=(double*)&this->y_ptr[this->element_number-2];
        }

        if((x_point<this->x_ptr[1]) && (x_point>this->x_ptr[this->element_number-2]))
        {
            for(int i=1;i<this->element_number-2;i++)
            {
                if((x_point<=this->x_ptr[i]) && (x_point>=this->x_ptr[i+1]))
                {
                    param.x=(double*)&this->x_ptr[i];
                    param.y=(double*)&this->y_ptr[i];
                }
            }
        }

    }


     mpfit(mpfit_routine,2,  this->routine_coef_numb, param_arr, 0, &config, (void *) &param, &result);

    return    routine(x_point,param_arr);

}

int mpfit_routine(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
    int i;
    struct vars_struct *v = static_cast<struct vars_struct*>(vars);
    double *x, *y, *ey, f;

    x = v->x;
    y = v->y;
    ey = v->ey;

    for (i=0; i<m; i++) {
        f= v->func(x[i],p);
        dy[i] = (y[i] - f)/ey[i];
    }

    return 0;
}
