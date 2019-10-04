/*
 *  Macros d'aide au debug.
 *  Définir la macro Use_Spy pour utiliser les fonctionnalités.
 *  Définir la macro SpyLevel fournissant le niveau des affichages.
 *
*/

#include <iostream>
#include <string>

#ifndef _Spy_Spy_H
#define _Spy_Spy_H

#ifndef SpyLevel
#define SpyLevel 10
#endif

#ifdef Use_Spy

#define SpyBlock(value) if(value<SpyLevel)
#define Spy_ConditionnalCode(code) code

#else

#define SpyBlock(value) if(false)
#define Spy_ConditionnalCode(code)

#endif

namespace Spy
{
	inline void verify(bool value, std::string const & message)
	{
		Spy_ConditionnalCode( if(!value) { std::cout<<message<<std::endl ; } ) 
	}
} 

#undef Spy_ConditionnalCode

#endif
