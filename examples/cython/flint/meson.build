project('python-flint-direct-example', 'cython', 'c', meson_version: '>=1.3')
py = import('python').find_installation(pure: false)
flint_dep = dependency('flint')
flint_include = run_command(
  py, '-c', 'import flint; print(flint.get_include())', check: true,
).stdout().strip()
py.extension_module('_flint_direct_example', '_flint_direct_example.pyx',
                    dependencies: flint_dep,
                    include_directories: include_directories(flint_include),
                    install: true)
test_env = environment()
test_env.prepend('PYTHONPATH', meson.current_build_dir())
test('direct FLINT example', py, args: files('test_flint.py'),
     env: test_env)
