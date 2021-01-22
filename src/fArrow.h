#ifndef __Arrow_fsh__
#define __Arrow_fsh__

inline std::string Arrow_fsh() 
{
	std::string str = 
	"varying vec4 fColor;\n"\
	"\n"\
	"uniform sampler2D pic_tex;\n"\
	"\n"\
	"void main()\n"\
	"{\n"\
	"	gl_FragColor = fColor;\n"\
	"\n"\
	"\n"\
	"\n"\
	"}\n";
	return str;
}


#endif
