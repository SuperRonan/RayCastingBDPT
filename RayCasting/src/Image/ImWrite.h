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
#include <Image/MultiSample.h>

namespace Image
{
	class ImWrite : public ImageIO
	{

	protected:

		template <bool MAJOR>
		static bool writeEXR(std::string const& path, Image<Geometry::RGBColor, MAJOR> const& image)
		{
			float* channels[3] = { new float[image.size()], new float[image.size()], new float[image.size()] };
			std::vector<std::string> cnames = { "R", "G", "B" };
			OMP_PARALLEL_FOR
				for (long i = 0; i < image.height(); ++i)
				{
					for (long j = 0; j < image.width(); ++j)
					{
						for (int c = 0; c < 3; ++c)
						{
							channels[c][i * image.width() + j] = image(j, i)[c];
						}
					}
				}
			Imf::Header header(image.width(), image.height());

			for (int c = 0; c < 3; ++c)
				header.channels().insert(cnames[c], Imf::Channel(Imf::FLOAT));

			Imf::OutputFile file(path.c_str(), header);


			Imf::FrameBuffer frame_buffer;
			for (int c = 0; c < 3; ++c)
				frame_buffer.insert(cnames[c], Imf::Slice(Imf::FLOAT, (char*)channels[c], sizeof(float), sizeof(float) * image.width()));

			file.setFrameBuffer(frame_buffer);

			file.writePixels(image.height());



			for (int j = 0; j < 3; ++j)
				delete[] channels[j];
			return true;
		}



		static bool writePPM_full(std::string const& path, const uint8_t* buffer, size_t width, size_t height)
		{
			try
			{
				std::stringstream file;
				std::cout << "saving the image: " << path << std::endl;
				file << "P3\n";
				file << width << " " << height << "\n";
				file << "255\n";
				for (size_t i = 0; i < width; ++i)
				{
					for (size_t j = 0; j < height; ++j)
					{
						size_t index = j * width + i;
						short r = buffer[index * 3];
						short g = buffer[index * 3 + 1];
						short b = buffer[index * 3 + 2];
						file << r << " " << g << " " << b << "\n";
					}
				}

				std::ofstream the_file(path);
				the_file << file.str();
				the_file.close();
				return true;
			}
			catch (std::exception const& e)
			{
				std::cerr << "Exception: " << std::endl;
				std::cerr << e.what() << std::endl;
				return false;
			}
		}


		////////////////////////
		//writes a ppm of the buffer data in path
		////////////////////////
		static bool writePPM_buffer(std::string const& path, const unsigned char* data, size_t width, size_t height)
		{
			try
			{
				std::ofstream file(path);
				file << "P6\n";
				file << width << " " << height << " ";
				file << "255 ";
				file.write((const char*)data, width * height * 3);
				file.close();
				return true;
			}
			catch (std::exception const& e)
			{
				std::cerr << "Exception: " << std::endl;
				std::cerr << e.what() << std::endl;
				return false;
			}

		}
		template <bool MAJOR>
		static bool writePPM(std::string const& path, Image<Geometry::RGBColor, MAJOR> img)
		{
			unsigned char* data = new unsigned char[img.size() * 3];

			OMP_DYNAMIC_FOR
				for (long i = 0; i < img.size(); ++i)
				{
					Geometry::RGBColor color = img.m_data[i];
#ifdef MULT_10
					color = color * 10;
#endif
					unsigned char r = (unsigned char)(color[0] / (color[0] + 1) * 255);
					unsigned char g = (unsigned char)(color[1] / (color[1] + 1) * 255);
					unsigned char b = (unsigned char)(color[2] / (color[2] + 1) * 255);
					data[i * 3 + 0] = r;
					data[i * 3 + 1] = g;
					data[i * 3 + 2] = b;
				}

			bool res = writePPM_full(path, data, img.width(), img.height());
			delete[] data;
			return res;
		}






	protected:
		static std::string getPath()
		{
			std::string path = RESULT_FOLDER;
			std::cout << "Where do you want to save the result:" << std::endl;
			std::cout << path;
			std::string name;
			std::cin >> name;

			if (name == "0")
			{
				return "";
			}

			return path + name;
		}

	public:

		template <bool MAJOR>
		static bool write(Image<Geometry::RGBColor, MAJOR> const& img, std::string path)
		{
			bool res;
			std::cout << "Writing " << path << std::endl;
			if (extension(path) == ".ppm")
			{
				res = writePPM(path, img);
			}
			else if (extension(path) == ".exr")
			{
				res = writeEXR(path, img);
			}
			else
			{
				std::cerr << "Unsuported format " << extension(path) << std::endl;
				res = false;
			}
			if (res)
				std::cout << "Done!" << std::endl;
			else
				std::cout << "Fail!" << std::endl;
			return res;
		}


		template <bool MAJOR>
		static bool write(Image<Geometry::RGBColor, MAJOR> const& img)
		{
			std::string path = getPath();
			if (path.empty())
			{
				std::cout << "aborting!" << std::endl;
				return false;
			}
			bool res = write(img, path);
			return res;
		}



		template <bool MAJOR>
		static bool write(Image<MultiSample<Geometry::RGBColor>, MAJOR> const& img)
		{
			Image < Geometry::RGBColor> _img(img.width(), img.height());
			OMP_DYNAMIC_FOR
				for (long i = 0; i < img.size(); ++i)
				{
					Geometry::RGBColor color = img.m_data[i].mean();
					_img.m_data[i] = color;
				}

			bool res = write(_img);
			return res;
		}

		template <bool MAJOR>
		static bool write(Image<Geometry::RGBColor, MAJOR> const& img, double f)
		{
			Image < Geometry::RGBColor> _img(img.width(), img.height());
			OMP_DYNAMIC_FOR
				for (long i = 0; i < img.size(); ++i)
				{
					Geometry::RGBColor color = img.m_data[i];
					_img.m_data[i] = color * f;
				}

			bool res = write(_img);
			return res;
		}


	};
}