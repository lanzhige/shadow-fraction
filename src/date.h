#ifndef CALCULATEMRT_DATE_H_
#define CALCULATEMRT_DATE_H_

#include <iostream>
#include <string>

namespace calculate {
class Date {
public:
	Date(const Date &);

	Date(const std::string &);

	Date& operator= (const Date& d);

	~Date() {}

	// if the year is a leap year
	bool isLeap(int year);

	// return the number of days in month.year
	int getDaysInMonth(int year, int month);

	bool operator== (const Date& d) const;
	bool operator!= (const Date& d) const;
	bool operator> (const Date& d);
	bool operator>= (const Date& d);
	bool operator< (const Date& d);
	bool operator<= (const Date& d);

	Date operator+ (int minute);
	Date& operator+= (int minute);
	Date operator- (int minute);
	Date& operator-= (int minute);
	Date operator++ (int);
	Date& operator++ ();
	Date operator-- (int);
	Date& operator-- ();

	/* get the minute between dates */
	int operator- (const Date& d);

	/* abandoned function */
	void nextHalfHour();

	/* return the sequence number of the day in the year */
	void getDayOfYear();


	std::string to_string();
  std::string to_utc();

	int dayOfYear;

	int _year;
	int _month;
	int _day;
	int _hour;
	int _minute;
	double sunHours;
	std::string _tzd;
	std::string str;

private:
	void update();
	void downdate();

	const int monthday[13] = { 31, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
};
}

#endif // !CALCULATEMRT_DATE_H_
