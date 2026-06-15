#include <stdio.h>

#include "flint/mpoly.h"
#include "flint/gr.h"
#include "flint/gr_mpoly.h"

#define TRACE(n) do { printf("TRACE %d\n", (n)); fflush(stdout); } while (0)

int main(int argc, char *argv[])
{
    gr_ctx_t ctx_fmpzi;

    printf("before...\n"); fflush(stdout);

    gr_ctx_init_fmpzi(ctx_fmpzi);

    printf("...after\n"); fflush(stdout);

    return 0;
}
