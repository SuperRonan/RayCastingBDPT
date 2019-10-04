#pragma once

#include <settings.h>
#include <string>
#include <Image/Image.h>
#include <Geometry/RGBColor.h>
#include <iostream>
#include <iosfwd>
#include <fstream>
#include <Image/ImageIO.h>
#include <ImfHeader.h>
#include <ImfOutputFile.h>
#include <ImfStringAttribute.h>
#include <ImfArray.h>
#include <ImfRgbaFile.h>
#include <ImfRgba.h>
#include <ImfChannelList.h>
#include <ImfInputFile.h>


namespace Image
{
	class ImRead: public ImageIO
	{
	protected:

	public:


		static bool readEXR(std::string const& path, Image<Geometry::RGBColor> & res)
		{
			Imf::InputFile file(path.c_str());
			Imath_2_3::Box2i dw = file.header().dataWindow();
			int width = dw.max.x - dw.min.x;
			int height = dw.max.y - dw.min.y;
			res.resize(width, height);
			float* channels[3];
			for (int c = 0; c < 3; ++c)
			{
				channels[c] = new float[res.size()];
			}
			std::string cnames[] = { "R", "G", "B" };

			Imf::FrameBuffer frame_buffer;

			for (int c = 0; c < 3; ++c)
			{
				frame_buffer.insert(cnames[c], Imf::Slice(
					Imf::FLOAT, (char*) (channels[c]), sizeof(float), sizeof(float) * width, 1, 1, 0.0)
				);

			}

			file.setFrameBuffer(frame_buffer);
			file.readPixels(dw.min.y, dw.max.y);

			for (size_t i = 0; i < res.width(); ++i)
			{
				for (size_t j = 0; j < res.height(); ++j)
				{
					for (int c = 0; c < 3; ++c)
					{
						res[i][j][c] = channels[c][j * width + i];
					}
				}
			}
			return true;
		}




		bool read(std::string path, Image<Geometry::RGBColor>& res)
		{
			std::string ext = extension(path);
			if (ext == ".exr")
			{

			}
			else if (ext == ".ppm")
			{

			}
			else if(ext == ".png")
			{

			}
			else
			{
				std::cerr << "Unknown file type: " << path << std::endl;
				return false;
			}
		}
	};
}