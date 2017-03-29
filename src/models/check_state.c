#if !defined(_MSC_VER)
#include <config.h>
#else 
#include <config_msvc.h>
#endif

#include <assert.h>
#include <math.h>

#include <structs.h>
#include <models/check_state.h>


//Type 40 / 261 / 262
int dam_check_qvs(
    double *y, unsigned int num_dof,
    const double * const params, unsigned int num_params,
    const double * const global_params, unsigned int num_global_params,
    QVSData *qvs,
    bool has_dam,
    void *user)
{
    unsigned int i, iterations;
    double S = y[1];

    if (!has_dam)
        return -1;

    iterations = qvs->n_values - 1;
    for (i = 0; i<iterations; i++)
    {
        if (qvs->points[i][0] <= S && S < qvs->points[i + 1][0])
            return i;
    }

    return i;
}