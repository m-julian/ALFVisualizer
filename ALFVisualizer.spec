# -*- mode: python ; coding: utf-8 -*-

block_cipher = None


a = Analysis(['ALFVisualizer.py'],
             pathex=['D:\\VMShared\\pyvista_alf'],
             binaries=[],
             datas=[("ALFVisualizer.ui", ".")],
             hiddenimports=['vtkmodules', 'vtkmodules.all', 'vtkmodules.util.numpy_support', 'vtkmodules.numpy_interface', 'vtkmodules.numpy_interface.dataset_adapter','vtkmodules.qt', 'vttmodules.util','vttmodules.vtkCommonCore','vttmodules.vtkCommonKitPython','vtkmodules.qt.QVTKRenderWindowInteractor'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='ALFVisualizer',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=False )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               upx_exclude=[],
               name='ALFVisualizer')
