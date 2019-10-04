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
		size_t m_width ;
		/// \brief	Window height.
		size_t m_height ;

		int m_scale;

		size_t m_ext_w;

		size_t m_ext_h;

		
		/// \brief	The window.
		SDL_Window * m_window;
		/// \brief	The renderer.
		SDL_Renderer * m_renderer;

		SDL_Surface * m_frame_buffer;

		SDL_Texture * m_tex_buffer;

		

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
		Visualizer(int width, int height, int scale = 1):
			m_width(width), 
			m_height(height), 
			m_scale(scale),
			m_ext_w(width/scale),
			m_ext_h(height/scale),
			m_frame_buffer(NULL)
		{
			Uint32 rmask, gmask, bmask, amask;

			/* SDL interprets each pixel as a 32-bit number, so our masks must depend
			   on the endianness (byte order) of the machine */
#if SDL_BYTEORDER == SDL_BIG_ENDIAN
			rmask = 0xff000000;
			gmask = 0x00ff0000;
			bmask = 0x0000ff00;
			amask = 0x000000ff;
#else
			rmask = 0x000000ff;
			gmask = 0x0000ff00;
			bmask = 0x00ff0000;
			amask = 0xff000000;
#endif
			if (SDL_Init(SDL_INIT_VIDEO) < 0)
			{
				::std::cerr << "Critical error" << ::std::endl;
				::std::cerr << "SDL_Init problem: " << SDL_GetError() << ::std::endl;
				exit(1);
			}
			atexit(SDL_Quit);

			m_window = SDL_CreateWindow("RTxRT", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, m_width, m_height, 0);

			m_renderer = SDL_GetRenderer(m_window);
			
			
			//SDL_CreateWindowAndRenderer(width * scale, height * scale, 0, &m_window, &m_renderer);
			
							
			//m_frame_buffer = SDL_CreateRGBSurface(0, m_width * m_scale, m_height * m_scale, 32, rmask, gmask, bmask, amask);
			m_frame_buffer = SDL_GetWindowSurface(m_window);
			if (m_frame_buffer == NULL)
			{
				std::cerr << "Critical Error!" << std::endl;
				std::cerr << "Could not create the frame buffer!" << std::endl;
				std::cerr << SDL_GetError() << std::endl;
				exit(1);
			}

			/*
			m_tex_buffer = SDL_CreateTextureFromSurface(m_renderer, m_frame_buffer);
			if (m_tex_buffer == NULL)
			{
				std::cerr << "Error, could not create the texture!" << std::endl;
				std::cerr << SDL_GetError() << std::endl;
				exit(1);
			}
			*/
		}

		~Visualizer()
		{
			SDL_FreeSurface(m_frame_buffer);
			SDL_DestroyRenderer(m_renderer);
			SDL_DestroyWindow(m_window);
			SDL_DestroyTexture(m_tex_buffer);
		}



		inline void update_scale(int new_scale)
		{
			if (new_scale != m_scale)
			{
				m_scale = new_scale;
				m_ext_w = m_width / m_scale;
				m_ext_h = m_height / m_scale;
			}
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
		size_t width() const
		{ return m_ext_w ; }

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
		size_t height() const
		{ return m_ext_h ; }


		size_t frameBufferWidth()const
		{
			return m_width;
		}

		size_t frameBufferHeight()const
		{
			return m_height;
		}



		inline Uint32 make_color(unsigned char r, unsigned char g, unsigned char b)const
		{
			return SDL_MapRGB(m_frame_buffer->format, r, g, b);
		}


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
		inline void plot(int x, int y, unsigned char r, unsigned char g, unsigned char b) 
		{
			//SDL_SetRenderDrawColor(m_renderer, r, g, b, SDL_ALPHA_OPAQUE);
			SDL_Rect rect;

			rect.x = x * m_scale;
			rect.y = y * m_scale;
			rect.w = m_scale;
			rect.h = m_scale;

			SDL_FillRect(m_frame_buffer, &rect, make_color(r, g, b));
		}

#ifdef MULT_10
		inline void plot(int x, int y, Geometry::RGBColor color){
			color = color * 10;
#else
		inline void plot(int x, int y, const Geometry::RGBColor& color){
#endif
			// A Simple tone mapper 
			unsigned char r = (unsigned char)(color[0] / (color[0] + 1) * 255);
			unsigned char g = (unsigned char)(color[1] / (color[1] + 1) * 255);
			unsigned char b = (unsigned char)(color[2] / (color[2] + 1) * 255);
			// Draws the pixel
			plot(x, y, r, g, b);
		}



		/*
		void draw_line_renderer(Math::Vector2f const& a, Math::Vector2f const& b, RGBColor const& color)
		{
			
		}
		*/

		enum KeyboardRequest {none, done, save};

		KeyboardRequest update_keyboard()
		{
			events.clear();
			KeyboardRequest res = none;
			SDL_Event event;
			while (SDL_PollEvent(&event)) {
				switch (event.type) {
				case SDL_KEYDOWN:
					if (event.key.keysym.sym == SDLK_ESCAPE)
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
					if (event.key.keysym.sym == SDLK_RETURN)
					{
						res = done;
					}
					else if (event.key.keysym.sym == SDLK_o)
					{
						res = save;
					}
					else
					{
						events.push_back(event);
					}
					break;
				case SDL_QUIT:
					exit(0);
					break;
				default:
					break;
				}
			}/*while*/

			return res;
		}


		bool update_renderer()
		{
			SDL_RenderPresent(m_renderer);
			return update_keyboard();
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
		KeyboardRequest update()
		{
			if (visible())
			{
				SDL_UpdateWindowSurface(m_window);
			}
			return update_keyboard();
		}


		template <class out_t>
		void print_info(out_t & out)const
		{
			out << "Visualizer: \n";
			out << "window resolution: " << m_width << "x" << m_height<<"\n";
			out << "Computing resolution: " << m_ext_w << "x" << m_ext_h<<"\n";
			out << "Scale: " << m_scale << std::endl;
		}



		void clean(unsigned char r=0, unsigned char g=0, unsigned char b=0)
		{
			SDL_Rect rect;

			rect.x = 0;
			rect.y = 0;
			rect.w = m_width;
			rect.h = m_height;
			SDL_FillRect(m_frame_buffer, &rect, make_color(r, g, b));
		}

		/////////////////////////////////
		//Warning ineficient function, use for debug
		/////////////////////////////////
		void fill(Geometry::RGBColor const& color)
		{
			for (size_t j = 0; j < height(); ++j)
			{
				for (size_t i = 0; i < width(); ++i)
				{
					plot(i, j, color);
				}
			}
		}

		const unsigned char* frameBufferData()const
		{
			return (const unsigned char*) m_frame_buffer->pixels;
		}

		bool visible()const
		{
			return !(SDL_GetWindowFlags(m_window) & SDL_WINDOW_MINIMIZED);
		}

	} ;
}

#endif
