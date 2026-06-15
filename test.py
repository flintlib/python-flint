from flint.types._gr import gr_fmpzi_ctx, gr_gr_mpoly_ctx
ctx = gr_gr_mpoly_ctx.new(gr_fmpzi_ctx, ["x", "y"])
print(repr(ctx.gens()))
x, y = ctx.gens()
p = (x + y)**2
print(repr(p))
