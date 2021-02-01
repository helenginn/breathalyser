project('splitseq', 'cpp', 'c')
qt5 = import('qt5')

qt5_dep = dependency('qt5', modules: ['Core', 'Gui', 'Widgets'], required : true)
dep_gl = dependency('gl', required : true)
png_dep = dependency('libpng')
boost_dep = dependency('boost')
fftw_dep = dependency('fftw3f')
dep_vag = dependency('libvag')
dep_vgeom = dependency('libvgeom', required: false)
dep_ccp4 = dependency('libccp4')
dep_c4x = dependency('libcluster4x')
dep_ica = dependency('libica')
helen3d_dep = dependency('helen3d')
helencore_dep = dependency('helencore')

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)

# Run the version and hash generator
hash_script = find_program('get_hash.sh')

gen_src = custom_target('commit-output',
                        input : ['src/commit.h.pre'],
                        output : ['commit.h'],
                        command : [hash_script, '@INPUT@', '@OUTPUT@'],
                        build_by_default: true,
                        build_always_stale: true)

# Graphical user interface
# (including all the Qt and OpenGL stuff)

# Pre-processing by Qt5

moc_files = qt5.preprocess(moc_headers : [
'src/CoupleDisplay.h', 
'src/Difference.h', 
'src/FastaMaster.h', 
'src/FastaGroup.h', 
'src/Fasta.h', 
'src/Main.h', 
'src/MutationWindow.h', 
'src/LoadFastas.h', 
'src/LoadStructure.h', 
'src/SlidingWindow.h', 
'src/SequenceView.h', 
'src/Segment.h', 
],
		moc_extra_arguments: ['-DMAKES_MY_MOC_HEADER_COMPILE'])

executable('splitseq', gen_src, moc_files,
'src/Arrow.cpp', 
'src/CoupleDisplay.cpp', 
'src/DiffDisplay.cpp', 
'src/Difference.cpp', 
'src/Ensemble.cpp', 
'src/Fasta.cpp', 
'src/FastaGroup.cpp', 
'src/FastaMaster.cpp', 
'src/LoadFastas.cpp', 
'src/LoadStructure.cpp', 
'src/Main.cpp', 
'src/MutationWindow.cpp', 
'src/MyDictator.cpp', 
'src/Segment.cpp', 
'src/SequenceView.cpp', 
'src/SlidingWindow.cpp', 
'src/StructureView.cpp', 
'src/WidgetFasta.cpp', 
'src/main.cpp', 
cpp_args: ['-std=c++11'], 
dependencies : [ helencore_dep, helen3d_dep, qt5_dep, dep_gl, png_dep, dep_vag, dep_vgeom, dep_ccp4, dep_c4x, boost_dep ], install: true)

