#pragma once

#include <iostream>
#include <windows.h>
#include <utils.h>
#include <iomanip>

class ProgressReporter
{
protected:

	size_t m_total;

	int m_width;

	LARGE_INTEGER m_frequency;
	LARGE_INTEGER m_begin, m_current;

	double m_time=0;

	bool m_new_line = false;

public:

	ProgressReporter(bool new_line=true):
		m_new_line(new_line)
	{}

	void start(size_t total)
	{
		m_total = total;
		CONSOLE_SCREEN_BUFFER_INFO console_info;
		GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &console_info);
		m_width = console_info.srWindow.Right - console_info.srWindow.Left + 1;
		
		QueryPerformanceFrequency(&m_frequency);
		QueryPerformanceCounter(&m_begin);
		m_time = 0;
	}


	void report(size_t n, int display_offset=0)
	{
		std::cout << std::setprecision(3) << std::fixed;
		QueryPerformanceCounter(&m_current);
		double current_time = (double)(m_current.QuadPart - m_begin.QuadPart) / (double)m_frequency.QuadPart;
		double tick_time = current_time - m_time;
		double remaining_time = (current_time / n) * (m_total - n);
		if (!m_new_line)
			std::cout << '\r';
		std::cout << n+display_offset << " / " << m_total << ": " 
			<< current_time << "s + " << remaining_time << "s = " 
			<< current_time + remaining_time << "s ("<< tick_time << "s), "<<percent(n, m_total) << "            ";
		if (m_new_line)
			std::cout << '\n';
		m_time = current_time;
	}

	double time()const
	{
		return m_time;
	}

	double finish()
	{
		QueryPerformanceCounter(&m_current);
		double current_time = (double)(m_current.QuadPart - m_begin.QuadPart) / (double)m_frequency.QuadPart;
		std::cout << "\Done! Total: " << current_time << "s, avg: " << current_time / m_total << "s" << std::endl;
		m_time = current_time;
		return current_time;
	}
};
