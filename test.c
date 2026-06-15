#include <stdio.h>

#include "flint/mpoly.h"
#include "flint/gr.h"
#include "flint/gr_mpoly.h"

#define TRACE(n) do { printf("TRACE %d\n", (n)); fflush(stdout); } while (0)

int main(int argc, char *argv[])
{
    int err;
    const char *names[] = {"x", "y"};
    char *pretty;
    gr_ptr x, y, r1, r2;

    gr_ctx_t ctx_fmpzi, ctx_mpoly;
    TRACE(1);
    gr_ctx_init_fmpzi(ctx_fmpzi);
    TRACE(2);
    gr_ctx_init_gr_mpoly(ctx_mpoly, ctx_fmpzi, 2, ORD_LEX);
    TRACE(3);
    err = gr_ctx_set_gen_names(ctx_mpoly, names);
    TRACE(4);
    if(err != GR_SUCCESS) {
        TRACE(5);
        printf("Cannot set generator names\n");
        fflush(stdout);
        TRACE(6);
        return 1;
    }
    TRACE(7);

    GR_TMP_INIT4(x, y, r1, r2, ctx_mpoly);
    TRACE(8);

    gr_vec_t gens;
    TRACE(9);
    gr_vec_init(gens, 0, ctx_mpoly);
    TRACE(10);
    err = gr_gens(gens, ctx_mpoly);
    TRACE(11);
    if(err != GR_SUCCESS) {
        TRACE(12);
        printf("Cannot get generators\n");
        fflush(stdout);
        TRACE(13);
        return 1;
    }
    TRACE(14);
    err |= gr_set(x, gr_vec_entry_ptr(gens, 0, ctx_mpoly), ctx_mpoly);
    TRACE(15);
    err |= gr_set(y, gr_vec_entry_ptr(gens, 1, ctx_mpoly), ctx_mpoly);
    TRACE(16);
    gr_vec_clear(gens, ctx_mpoly);
    TRACE(17);

    err |= gr_add(r1, x, y, ctx_mpoly);
    TRACE(18);
    err |= gr_pow_si(r2, r1, 2, ctx_mpoly);
    TRACE(19);

    if(err != GR_SUCCESS) {
        TRACE(20);
        printf("Arithmetic failed\n");
        fflush(stdout);
        TRACE(21);
        return 1;
    }
    TRACE(22);

    err = gr_get_str(&pretty, r2, ctx_mpoly);
    TRACE(23);

    if(err != GR_SUCCESS) {
        TRACE(24);
        printf("Pretty printing failed\n");
        fflush(stdout);
        TRACE(25);
        return 0;
    }
    TRACE(26);

    printf("(x + y)^2 = %s\n", pretty);
    fflush(stdout);
    TRACE(27);

    return 0;
}
