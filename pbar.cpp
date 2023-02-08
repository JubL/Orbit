#include "pbar.h"

Pbar::Pbar()
{
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &_winsize);
//	_lastcol = 1;
}


void Pbar::pBarInit()
{
	_lastCol = 1;
	for(int32_t i = 0; i < _winsize.ws_col; i++)
	{
		mvaddch(_winsize.ws_row - 1, i, ' ');
	}
	mvaddch(_winsize.ws_row - 1, 0, '[');
        mvaddch(_winsize.ws_row - 1, _winsize.ws_col - 10, ']');
	refresh();
	curs_set(0);
}

void Pbar::progressBar(double i)
{
	int32_t currentCol = std::lround((double)(_winsize.ws_col - 11) * i / 100.0 + 1);
	for(; _lastCol < currentCol - 1; _lastCol++)
	{
		mvaddch(_winsize.ws_row - 1, _lastCol, '=');
	}
	mvaddch(_winsize.ws_row - 1, _lastCol, i > 99.99 ? '=' : '>');
	std::ostringstream o;
	o << std::setw(6) << std::setfill(' ') << std::fixed << std::setprecision(2) << i << '%';
	std::string percent = o.str();
	mvaddch(_winsize.ws_row - 1, _winsize.ws_col - 8, percent.at(0));
	mvaddch(_winsize.ws_row - 1, _winsize.ws_col - 7, percent.at(1));
	mvaddch(_winsize.ws_row - 1, _winsize.ws_col - 6, percent.at(2));
	mvaddch(_winsize.ws_row - 1, _winsize.ws_col - 5, percent.at(3));
	mvaddch(_winsize.ws_row - 1, _winsize.ws_col - 4, percent.at(4));
	mvaddch(_winsize.ws_row - 1, _winsize.ws_col - 3, percent.at(5));
	mvaddch(_winsize.ws_row - 1, _winsize.ws_col - 2, percent.at(6));
	refresh();
}

void Pbar::showEstimatedTimeTillFinish(time_t& start, double advance, double settingsPosRow)
{
	time_t now;
	time_t estimatedEnd;
	time(&now);
	double diffTime = difftime(now, start);
	diffTime /= advance;
	estimatedEnd = start + (time_t)diffTime;
	char mbstr[100];
	std::strftime(mbstr, sizeof(mbstr), "%A %d. %B %Y %R %Z", std::localtime(&estimatedEnd));
	mvprintw(settingsPosRow + 21, 10, "%s", "Estimated time of finish:");
	mvprintw(settingsPosRow + 21, 36, "%s   ", mbstr);
}

void Pbar::printInfo(double atmTemp, double soilTemp, Star& hoststar, Planet& exoplanet, std::string filename, time_t& start, double advance)
{
	double settingsPosRow = _winsize.ws_row/2;
	double settingsPosColStar = 10;
	double settingsPosColExo = 55 + settingsPosColStar;

	if(advance && (int)(advance*10000) % 100 == 0) showEstimatedTimeTillFinish(start, advance, settingsPosRow); // advance == 0 => if(false)

	mvprintw(settingsPosRow + 9, 10, "%s", "Surface air temperature:");
	mvprintw(settingsPosRow + 9, 33, "%7.2f K", atmTemp);
	mvprintw(settingsPosRow + 10, 33, "%7.2f °C", atmTemp-273.15);

	mvprintw(settingsPosRow + 12, 10, "%s", "Surface temperature:");
	mvprintw(settingsPosRow + 12, 33, "%7.2f K", soilTemp);
	mvprintw(settingsPosRow + 13, 33, "%7.2f °C", soilTemp-273.15);

	mvprintw(settingsPosRow + 15, 10, "%s", "Output filename:");
	mvprintw(settingsPosRow + 15, 33, "%s", filename.c_str());

	mvprintw(settingsPosRow + 17, 10, "%s", "Current orbital phase:");
	mvprintw(settingsPosRow + 17, 33, "%.2f Pi", exoplanet.getPhaseAngle()/M_PI);

	mvprintw(settingsPosRow, settingsPosColStar, "%s", "Star:");
	mvprintw(settingsPosRow + 1, settingsPosColStar + 5, "%s", "Mass:");
	mvprintw(settingsPosRow + 1, settingsPosColStar + 20, "%.2f solar mass", hoststar.getMass(1));
	mvprintw(settingsPosRow + 2, settingsPosColStar + 5, "%s", "Radius:");
	mvprintw(settingsPosRow + 2, settingsPosColStar + 20, "%.2f solar radii", hoststar.getRadius(1));
	mvprintw(settingsPosRow + 3, settingsPosColStar + 5, "%s", "Temperature:");
	mvprintw(settingsPosRow + 3, settingsPosColStar + 20, "%.2f Kelvin", hoststar.getTemperature());

	mvprintw(settingsPosRow, settingsPosColExo, "%s", "Exoplanet:");
	mvprintw(settingsPosRow + 1, settingsPosColExo + 5, "%s", "Mass:");
	mvprintw(settingsPosRow + 1, settingsPosColExo + 24, "%.3f earth mass", exoplanet.getMass(1));
	mvprintw(settingsPosRow + 2, settingsPosColExo + 5, "%s", "Radius:");
	mvprintw(settingsPosRow + 2, settingsPosColExo + 24, "%.3f earth radii", exoplanet.getRadius(1));
	mvprintw(settingsPosRow + 3, settingsPosColExo + 5, "%s", "Albedo:");
	mvprintw(settingsPosRow + 3, settingsPosColExo + 24, "%.3f", exoplanet.getAlbedo());
	mvprintw(settingsPosRow + 4, settingsPosColExo + 5, "%s", "Eccentricity:");
	mvprintw(settingsPosRow + 4, settingsPosColExo + 24, "%.3f", exoplanet.getEccentricity());
	mvprintw(settingsPosRow + 5, settingsPosColExo + 5, "%s", "Semi major axis:");
	mvprintw(settingsPosRow + 5, settingsPosColExo + 24, "%.3f astronomical units", exoplanet.getSemiMajorAxis(1));
	mvprintw(settingsPosRow + 6, settingsPosColExo + 5, "%s", "Orbital period:");
	mvprintw(settingsPosRow + 6, settingsPosColExo + 24, "%.3f earth days", exoplanet.getOrbitalPeriod(1));
	mvprintw(settingsPosRow + 7, settingsPosColExo + 5, "%s", "H2O:");
	mvprintw(settingsPosRow + 7, settingsPosColExo + 24, "%.3e", exoplanet.getVolumeRatioH2O());
	mvprintw(settingsPosRow + 8, settingsPosColExo + 5, "%s", "CO2:");
	mvprintw(settingsPosRow + 8, settingsPosColExo + 24, "%.3e", exoplanet.getVolumeRatioCO2());
	mvprintw(settingsPosRow + 9, settingsPosColExo + 5, "%s", "CH4:");
	mvprintw(settingsPosRow + 9, settingsPosColExo + 24, "%.3e", exoplanet.getVolumeRatioCH4());
}

void Pbar::pBarStop()
{
	endwin();
}
