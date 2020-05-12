#pragma once

#include <Image/Image.h>
#include <string>

namespace Image
{
	class ImageIO
	{
	protected:

	public:


		static std::string extension(std::string const& file)
		{
			std::string res;
			for (int i = file.size() - 1; i >= 0; --i)
			{
				if (file[i] == '/' || file[i] == '\\')
					break;
				res += file[i];
				if (file[i] == '.')
				{
					std::reverse(res.begin(), res.end());
					return res;
				}
			}

			return "";
		}


	};
}