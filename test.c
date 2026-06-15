#include "flint/mpoly.h"
#include "flint/gr.h"
#include "flint/gr_mpoly.h"

int main(int argc, char *argv[])
{
    int err;
    const char *names[] = {"x", "y"};
    char *pretty;
    gr_ptr x, y, r1, r2;

    gr_ctx_t ctx_fmpzi, ctx_mpoly;
    gr_ctx_init_fmpzi(ctx_fmpzi);
    gr_ctx_init_gr_mpoly(ctx_mpoly, ctx_fmpzi, 2, ORD_LEX);
    err = gr_ctx_set_gen_names(ctx_mpoly, names);
    if(err != GR_SUCCESS) {
        printf("Cannot set generator names\n");
        return 1;
    }

    GR_TMP_INIT4(x, y, r1, r2, ctx_mpoly);

    gr_vec_t gens;
    gr_vec_init(gens, 0, ctx_mpoly);
    err = gr_gens(gens, ctx_mpoly);
    if(err != GR_SUCCESS) {
        printf("Cannot get generators\n");
        return 0;
    }
    err |= gr_set(x, gr_vec_entry_ptr(gens, 0, ctx_mpoly), ctx_mpoly);
    err |= gr_set(y, gr_vec_entry_ptr(gens, 1, ctx_mpoly), ctx_mpoly);
    gr_vec_clear(gens, ctx_mpoly);

    err |= gr_add(r1, x, y, ctx_mpoly);
    err |= gr_pow_si(r2, r1, 2, ctx_mpoly);

    if(err != GR_SUCCESS) {
        printf("Arithmetic failed\n");
        return 1;
    }

    err = gr_get_str(&pretty, r2, ctx_mpoly);

    if(err != GR_SUCCESS) {
        printf("Pretty printing failed\n");
        return 0;
    }

    printf("(x + y)^2 = %s\n", pretty);

    return 0;
}
