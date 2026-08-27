project('python-flint-capi-example', 'cython', 'c', meson_version: '>=1.3')
py = import('python').find_installation(pure: false)
flint_include = run_command(
  py, '-c', 'import flint; print(flint.get_include())', check: true,
).stdout().strip()
py.extension_module(
  '_flint_capi_example', '_flint_capi_example.pyx',
  include_directories: include_directories(flint_include),
  install: true,
)
test('capsule API example', py, args: files('test_capi.py'),
     env: {'PYTHONPATH': meson.current_build_dir()})
