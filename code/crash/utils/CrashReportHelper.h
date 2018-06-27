#include <sys/socket.h>
#ifdef __APPLE__
#include <net/if_dl.h>
#define AF_LOWLEVEL AF_LINK
#else
#include <linux/if_packet.h>
#define AF_LOWLEVEL AF_PACKET
#endif
#include <sys/types.h>
#include <ifaddrs.h>
#include <stdlib.h>
#include <string.h>
#include <net/if.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <set>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>

namespace casac {

    class CrashReportHelper {


        typedef std::vector<unsigned> tmac;

        public:
            std::string getUniqueId () {
            struct ifaddrs *ifap, *ifaptr;
            unsigned char *ptr;
            std::set<tmac> macs;

            if (getifaddrs(&ifap) == 0) {
                for(ifaptr = ifap; ifaptr != NULL; ifaptr = (ifaptr)->ifa_next) {

                    if ( ! (ifaptr->ifa_flags & IFF_LOOPBACK) &&
                        ! (ifaptr->ifa_flags & IFF_POINTOPOINT) ) { // don't count loopback
                        // awdl is an OSX specific interface and seems to get a new Mac address assigned at every reboot
                        // so we filter that out.
                        if (((ifaptr)->ifa_addr)->sa_family == AF_LOWLEVEL && strncmp("awdl",(ifaptr)->ifa_name,4)) {
                            #ifdef __APPLE__
                            ptr = (unsigned char *)LLADDR((struct sockaddr_dl *)(ifaptr)->ifa_addr);
                            #else
                            ptr = (unsigned char *)(((struct sockaddr_ll*)ifaptr->ifa_addr)->sll_addr);
                            #endif
                            macs.insert(tmac({*ptr, *(ptr+1), *(ptr+2), *(ptr+3), *(ptr+4), *(ptr+5)}));
                            /*
                            printf( "%s: %02x:%02x:%02x:%02x:%02x:%02x\n",
                                    (ifaptr)->ifa_name,
                                    *ptr, *(ptr+1), *(ptr+2), *(ptr+3), *(ptr+4), *(ptr+5) );
                            */
                        }
                    }
                }
                freeifaddrs(ifap);

                tmac comb = accumulate( macs.begin( ),macs.end( ),tmac({ 0,0,0,0,0,0 }),
                                        []( tmac acc, tmac mac ) {
                                            tmac result(acc.size());
                                            std::transform(acc.begin(), acc.end(), mac.begin(), result.begin(), std::plus<int>());
                                            return result;
                                        } );
                /*
                printf( "--------------------------------------------------\n");
                printf( "sum collected as a possible machine uid\n" );
                printf( "--------------------------------------------------\n");
                printf( "hash: %03x%03x%03x%03x%03x%03x\n",
                        comb[0],comb[1],comb[2],comb[3],comb[4],comb[5] );
                */
                std::stringstream sstream;
                sstream << std::hex << comb[0] << comb[1] << comb[2] << comb[3] << comb[4] << comb[5];
                std::string hash = sstream.str();
                return hash;
            } else {
                return "NA";
            }
        }
    };
}
