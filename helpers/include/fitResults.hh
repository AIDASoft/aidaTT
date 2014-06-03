#ifndef FITRESULTS_HH
#define FITRESULTS_HH

#include <ostream>

namespace aidaTT
{

    class fitResults
    {
        public:

            fitResults() : _valid(false)
            {
                ;
            };

            fitResults(bool valid, double chis, unsigned int ndf, double weightL, const trackParameters& tp)
                : _valid(valid), _chisquare(chis), _ndf(ndf), _lostweight(weightL), _estParams(tp) {};

            double chiSquare() const
            {
                return _chisquare;
            };
            unsigned int ndf() const
            {
                return _ndf;
            };
            double weightLost() const
            {
                return _lostweight;
            };
            bool areValid() const
            {
                return _valid;
            };

            trackParameters estimatedParameters() const
            {
                return _estParams;
            };

            void setResults(bool v, double chs, unsigned int n, double wl, const trackParameters& tp)
            {
                _valid = v;
                _chisquare = chs;
                _ndf = n;
                _lostweight = wl;
                _estParams = tp;

            };
        private:
            bool _valid;
            double _chisquare;
            unsigned int _ndf;
            double _lostweight;

            trackParameters _estParams;
    };

    inline std::ostream & operator << (std::ostream & os, const fitResults& fR)
    {
        os << " [fitResults]: { results are valid? " << fR.areValid() << "} , { chi^2/ndf : " << fR.chiSquare() << "/" << fR.ndf() <<  " } "
           << ", ( lost weight : " << fR.weightLost() << " with the track parameters: " << fR.estimatedParameters();
        return os ;
    }

}
#endif // FITRESULTS_HH
