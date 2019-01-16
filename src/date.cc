#include "date.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

namespace calculate {
	Date::Date(const Date& d)
		:_year(d._year)
		, _month(d._month)
		, _day(d._day)
		, _hour(d._hour)
		, _minute(d._minute)
		, _tzd(d._tzd)
		, dayOfYear(d.dayOfYear)
	{}

	Date::Date(const std::string &utc) {
		if (utc.length() == 22) {
			_year = std::stoi(utc.substr(0, 4));
			_month = std::stoi(utc.substr(5, 2));
			_day = std::stoi(utc.substr(8, 2));
			_hour = std::stoi(utc.substr(11, 2));
			_minute = std::stoi(utc.substr(14, 2));
			_tzd = utc.substr(16, 6);
		} else {
			std::cout << "unresolved date format: " << utc << std::endl;
		}
    if (_tzd[0] == '-') {
      _hour -= std::stoi(_tzd);
      _minute += std::stoi(_tzd.substr(4));
      update();
    } else {
      _hour -= std::stoi(_tzd);
      _minute -= std::stoi(_tzd.substr(4));
      downdate();
    }
		getDayOfYear();
	}

	Date& Date::operator= (const Date& d)
	{
		if (this != &d) {
			_year = d._year;
			_month = d._month;
			_day = d._day;
			_hour = d._hour;
			_minute = d._minute;
			_tzd = d._tzd;
			dayOfYear = d.dayOfYear;
		}
		return *this;
	}

	bool Date::isLeap(int year) {
		if (year % 4 == 0 && year % 100 != 0 || year % 400 == 0) {
			return true;
		}
		return false;
	}

	int Date::getDaysInMonth(int year, int month)	{
		if (month == 2) {
			if (isLeap(year)) return 29;
			else return 28;
		}
		else {
			return monthday[month];
		}
	}

	void Date::getDayOfYear() {
		int count = 0;
		for (int i = 0; i < _month; i++) {
			if (i == 2) {
				if (isLeap(_year)) {
					count += 29;
				}
				else count += 28;
			}
			else {
				count += monthday[i];
			}
		}
		count += _day;
		dayOfYear = count-1;
	}

	bool Date::operator== (const Date& d) const	{
		return _year == d._year
			&& _month == d._month
			&& _day == d._day
			&& _hour == d._hour
			&& _minute == d._minute
			&& _tzd == d._tzd;
	}

	bool Date::operator!= (const Date& d) const	{
		return !(*this == d);
	}
	bool Date::operator> (const Date& d) {
		if (_year > d._year) {
			return true;
		}
		if (_year == d._year) {
			if (_month > d._month) {
				return true;
			}
			if (_month = d._month) {
				if (_day > d._day) {
					return true;
				}
				if (_day == d._day) {
					if (_hour > d._hour) {
						return true;
					}
					if (_hour == d._hour) {
						if (_minute > d._minute) {
							return true;
						}
					}
				}
			}
		}
		return false;
	}
	bool Date::operator>= (const Date& d) {
		return (*this > d) || (*this == d);
	}

	bool Date::operator< (const Date& d) {
		return !(*this >= d);
	}

	bool Date::operator<= (const Date& d) {
		return !(*this > d);
	}

	Date Date::operator+ (int minute) {
		if (minute < 0) {
			minute = -minute;
			return (*this - minute);
		}
		Date temp(*this);
		if (minute == 0) {
			return temp;
		}
		temp._minute += minute;
		temp.update();
		return temp;
	}

	Date& Date::operator+= (int minute)
	{
		Date temp(*this);
		temp = (*this) + minute;
		*this = temp;
		return *this;
	}

	Date Date::operator- (int minute) {
		this->_minute -= minute;
		downdate();

		return *this;
	}

	Date& Date::operator-= (int minute)
	{
		Date temp(*this);
		temp = (*this) - minute;
		*this = temp;
		return *this;
	}
	Date Date::operator++ (int)
	{
		Date tmp(*this);
		(*this) += 1;
		return tmp;
	}
	Date& Date::operator++ ()
	{
		(*this) += 1;
		return *this;
	}
	Date Date::operator-- (int)
	{
		Date tmp(*this);
		(*this) -= 1;
		return tmp;
	}
	Date& Date::operator-- ()
	{
		(*this) -= 1;
		return *this;
	}

	/*
	* return the minutes between two dates
	*/

	int Date::operator- (const Date &d){
		int result = 0;
		result += this->_minute - d._minute;
		result += 60 * (this->_hour - d._hour);
		result += 24 * 60 * (this->_day - d._day);
		//jan.4th apr.1
		//TODO: year month smaller
		for (int i = d._month; i < this->_month; i++) {
			result += getDaysInMonth(d._year, i) * 24 * 60;
		}
		
		return result;
	}

	

	void Date::nextHalfHour() {
		_minute += 30;
		update();
	}

	void Date::update() {
		while (_minute >= 60) {
			_minute -= 60;
			_hour++;
		}
		while (_hour >= 24) {
			_hour -= 24;
			_day++;
		}
		while (_month > 12) {
			_month -= 12;
			_year++;
		}
		while (_day >= getDaysInMonth(_year, _month)) {
			_day -= getDaysInMonth(_year, _month);
			_month++;
			if (_month > 12) {
				_month -= 12;
				_year++;
			}
		}
	}

	void Date::downdate(){
		while (_minute < 0){
			_minute += 60;
			_hour--;
		}
		while (_hour<0){
			_hour += 24;
			_day--;
		}

		while (_day <= 0) {
			_month--;
			while (_month < 1) {
				_month += 12;
				_year -- ;
			}
			_day += getDaysInMonth(_year, _month);
		}
	}

	std::string Date::to_string() {
		std::stringstream ss;
		ss.fill('0');
		_tzd[3] = '-';
    Date temp(*this);
    if (_tzd[0] == '+') {
      _hour += std::stoi(_tzd);
      _minute += std::stoi(_tzd.substr(4));
      update();
    } else {
      _hour += std::stoi(_tzd);
      _minute -= std::stoi(_tzd.substr(4));
      downdate();
    }
		ss << std::setw(4) << _year << "-"
			<< std::setw(2) << _month << "-"
			<< std::setw(2) << _day << "T"
			<< std::setw(2) << _hour << "-"
			<< std::setw(2) << _minute
			<< _tzd;
		str = ss.str();
    *this = temp;
		return str;
	}

  std::string Date::to_utc() {
    std::stringstream ss;
    ss.fill('0');
    _tzd[3] = ':';
    Date temp(*this);
    if (_tzd[0] == '+') {
      _hour += std::stoi(_tzd);
      _minute += std::stoi(_tzd.substr(4));
      update();
    } else {
      _hour += std::stoi(_tzd);
      _minute -= std::stoi(_tzd.substr(4));
      downdate();
    }
    ss << std::setw(4) << _year << "-"
      << std::setw(2) << _month << "-"
      << std::setw(2) << _day << "T"
      << std::setw(2) << _hour << ":"
      << std::setw(2) << _minute
      << _tzd;
    str = ss.str();
    *this = temp;
    return str;
  }
}