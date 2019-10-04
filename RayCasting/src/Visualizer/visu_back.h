#ifndef _Visualizer_Visualizer_H
#define _Visualizer_Visualizer_H

#include <SDL.h>
//#include <SDL_draw.h>
#include <iostream>
#include <Geometry/RGBColor.h>

namespace Visualizer
{
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
		/// \brief	Windows width.
		int m_width ;
		/// \brief	Window height.
		int m_height ;
		/// \brief	The window.
		SDL_Window * m_window;
		/// \brief	The renderer.
		SDL_Renderer * m_renderer;

		

	public:
		std::vector<SDL_Event> events;
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
			if (SDL_Init(SDL_INIT_VIDEO) < 0)
			{
				::std::cerr << "Critical error" << ::std::endl;
				::std::cerr << "SDL_Init problem: " << SDL_GetError() << ::std::endl;
				exit(1);
			}
			atexit(SDL_Quit);

			SDL_CreateWindowAndRenderer(width, height, 0, &m_window, &m_renderer);
		}

		void resize(int w, int h)
		{
			SDL_SetWindowSize(m_window, w, h);
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
		void plot(int x, int y, unsigned char r, unsigned char g, unsigned char b, int scale=1) 
		{
			SDL_SetRenderDrawColor(m_renderer, r, g, b, SDL_ALPHA_OPAQUE);
			SDL_Rect rect;

			rect.x = x;
			rect.y = y;
			rect.w = scale;
			rect.h = scale;

			SDL_RenderFillRect(m_renderer, &rect);
			//SDL_RenderDrawPoint(m_renderer, x, y);
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
		void plot(int x, int y, const Geometry::RGBColor color, int scale=1)
		{
			// A Simple tone mapper 
			unsigned char r = (unsigned char)(color[0]/(color[0]+1)*255) ;
			unsigned char g = (unsigned char)(color[1]/(color[1]+1)*255) ;
			unsigned char b = (unsigned char)(color[2]/(color[2]+1)*255) ;
			// Draws the pixel
			plot(x,y,r,g,b, scale);
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
			events.clear();
			SDL_RenderPresent(m_renderer);
			SDL_Event event;
			while ( SDL_PollEvent(&event) ) {
				switch (event.type) {
				case SDL_KEYDOWN:
					if (event.key.keysym.sym == 27) // Escape key
					{
						::std::cout << "Do you really want to quit [y or Y = yes]? ";
						char answer;
						::std::cin >> answer;
						if (answer == 'y' || answer == 'Y')
						{
							exit(0);
						}
					}
					else
					{
						events.push_back(event);
					}
					break;
				case SDL_KEYUP:
					events.push_back(event);
					break;
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
