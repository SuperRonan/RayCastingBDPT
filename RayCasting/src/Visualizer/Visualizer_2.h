#ifndef _Visualizer_Visualizer_H
#define _Visualizer_Visualizer_H

#include <SDL.h>
#include <SDL_draw.h>
#include <iostream>
#include <Geometry/RGBColor.h>

namespace Visualizer
{
	/** \brief Classe permettant d'effectuer une rendu graphique en 2D */

	////////////////////////////////////////////////////////////////////////////////////////////////////
	/// \class	Visualizer
	///
	/// \brief	Opens a 2D rendering context.
	///
	/// \author	F. Lamarche, Université de Rennes 1
	/// \date	03/12/2013
	////////////////////////////////////////////////////////////////////////////////////////////////////
	class Visualizer
	{
	protected:
		/// \brief	The rendering context.
		SDL_Surface * screen;
		/// \brief	Windows width.
		int m_width ;
		/// \brief	Window height.
		int m_height ;

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Uint32 Visualizer::FastestFlags(Uint32 flags, unsigned int width, unsigned int height,
		/// 	unsigned int bpp)
		///
		/// \brief	Gets flags for SDL draw configuration.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/12/2013
		///
		/// \param	flags 	The flags.
		/// \param	width 	The width.
		/// \param	height	The height.
		/// \param	bpp   	The bits per pixel.
		///
		/// \return	The flags.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Uint32 FastestFlags(Uint32 flags, unsigned int width, unsigned int height, unsigned int bpp)
		{
			const SDL_VideoInfo *info;
			flags |= SDL_FULLSCREEN;
			info = SDL_GetVideoInfo();
			if ( info->blit_hw_CC && info->blit_fill ) 
			{
				flags |= SDL_HWSURFACE;
			}
			if ( (flags & SDL_HWSURFACE) == SDL_HWSURFACE ) 
			{
				if ( info->video_mem*1024 > (height*width*bpp/8) ) 
				{
					flags |= SDL_DOUBLEBUF;
				} 
				else 
				{
					flags &= ~SDL_HWSURFACE;
				}
			}

			return flags;
		}

	public:

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	Visualizer::Visualizer(int width, int height)
		///
		/// \brief	Constructor.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/12/2013
		///
		/// \param	width 	The width of the rendering window.
		/// \param	height	The height of the rendering window.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		Visualizer(int width, int height)
			: m_width(width), m_height(height) 
		{
			if(SDL_Init(SDL_INIT_VIDEO)<0) 
			{
				::std::cerr<<"Critical error"<<::std::endl ;
				::std::cerr<<"SDL_Init problem: "<<SDL_GetError()<<::std::endl;
				exit(1);
			}
			atexit(SDL_Quit);
			Uint32 videoflags=(SDL_SWSURFACE | SDL_ANYFORMAT) ;  //= FastestFlags(SDL_SWSURFACE | SDL_ANYFORMAT, width, height, 0) ;
			
			screen = SDL_SetVideoMode(width, height, 0, videoflags);

			if (!screen) 
			{
				::std::cerr<<"Critical error"<<::std::endl ;
				::std::cerr<<"I can not activate video mode "<<width<<"x"<<height<<" : "<<SDL_GetError();
				exit(2);
			}
			// Initializes SDL Draw
			Draw_Init();
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	int Visualizer::width() const
		///
		/// \brief	Gets the width of the rendering window.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		int width() const
		{ return m_width ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	int Visualizer::height() const
		///
		/// \brief	Gets the height of the rendering window.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/12/2013
		///
		/// \return	.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		int height() const
		{ return m_height ; }

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Visualizer::plot(int x, int y, unsigned char r, unsigned char g,
		/// 	unsigned char b) const
		///
		/// \brief	Change the color of a given pixel.
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/12/2013
		///
		/// \param	x	The x coordinate of the pixel.
		/// \param	y	The y coordinate of the pixel.
		/// \param	r	The red value [0..255].
		/// \param	g	The green value [0..255].
		/// \param	b	The blue value [0..255].
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void plot(int x, int y, unsigned char r, unsigned char g, unsigned char b) const
		{
			Uint32 color = SDL_MapRGB(screen->format, r, g, b) ;
			Draw_Pixel(screen, (Sint16)x, (Sint16)y, color) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Visualizer::plot(int x, int y, const RGBColor color) const
		///
		/// \brief	Plots with a simple tone mapper
		///
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	04/12/2013
		///
		/// \param	x	 	The x coordinate.
		/// \param	y	 	The y coordinate.
		/// \param	color	The color.
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void plot(int x, int y, const Geometry::RGBColor color) const
		{
			// A Simple tone mapper
			unsigned char r = color[0]/(color[0]+1)*255 ;
			unsigned char g = color[1]/(color[1]+1)*255 ;
			unsigned char b = color[2]/(color[2]+1)*255 ;
			// Maps the result into the rendering context
			Uint32 renderedColor = SDL_MapRGB(screen->format, r, g, b) ;
			Draw_Pixel(screen, (Sint16)x, (Sint16)y, renderedColor) ;
		}

		////////////////////////////////////////////////////////////////////////////////////////////////////
		/// \fn	void Visualizer::update()
		///
		/// \brief	Updates the rendering context.
		///
		/// \warning No modification is visible until this method is called. Application exists if any
		/// 		 key is pressed.
		/// 
		/// \author	F. Lamarche, Université de Rennes 1
		/// \date	03/12/2013
		////////////////////////////////////////////////////////////////////////////////////////////////////
		void update()
		{
			SDL_UpdateRect(screen, 0, 0, 0, 0);
			SDL_Event event;
			while ( SDL_PollEvent(&event) ) {
				switch (event.type) {
				case SDL_KEYDOWN:
					/*break;*/
				case SDL_QUIT:
					exit(0) ;
					break;
				default:
					break;
				}
			}/*while*/
		}
	} ;
}

#endif
