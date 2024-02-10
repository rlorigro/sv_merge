#include "Authenticator.hpp"
#include "Timer.hpp"
#include "bam.hpp"
#include "curl/curl.h"

using namespace sv_merge;

#include <functional>

using std::function;


static size_t WriteCallback(void *contents, size_t size, size_t nmemb, void *userp)
{
    ((std::string*)userp)->append((char*)contents, size * nmemb);
    return size * nmemb;
}


void curl_get(const string& url, string& result){
    result.clear();
    CURL *curl;
    CURLcode res;

    curl = curl_easy_init();
    if(curl) {
        curl_easy_setopt(curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, WriteCallback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &result);
        res = curl_easy_perform(curl);

        if(res != CURLE_OK) {
            fprintf(stderr, "curl_easy_perform() failed: %s\n",
                    curl_easy_strerror(res));
        }

        curl_easy_cleanup(curl);

        std::cout << result << std::endl;
    }
}


void parse_json_line(const string& line, pair<string,string>& result){
    result.first.clear();
    result.second.clear();
    size_t n_quotes = 0;

    for (auto c: line){
        if (c == '\"'){
            n_quotes++;
            continue;
        }

        if (n_quotes == 1){
            result.first += c;
        }

        if (n_quotes == 3){
            result.second += c;
        }
    }

    if (n_quotes != 4){
        throw runtime_error("ERROR: attempt to parse non-string value in JSON line: " + line);
    }
}


/**
 * only valid for simple STRING objects no nested objects, no ints, no floats...
 */
void for_item_in_json(const string& json, const function<void(const pair<string,string>& item)>& f){
    string token;
    pair<string,string> item;

    for (const auto c: json){
        if (c == ','){
            parse_json_line(token, item);
            f(item);
            token.clear();
        }
        else{
            token += c;
        }
    }
}


int main(){
    string command = "gcloud auth print-access-token";
    string token;

    // Will throw error if fails, otherwise returns token
    run_command(command, token);

    string url = "https://oauth2.googleapis.com/tokeninfo?access_token=" + token;
    string result;
    curl_get(url, result);

    for_item_in_json(result, [&](const pair<string,string>& item){
        if (item.first == "exp"){
            auto remaining_seconds = stoll(item.second);
        }
    });

    GoogleAuthenticator authenticator;
    string region_string = "chr1:10000000-10005000";
    path bam_path = "gs:/fc-b1dcc7b0-68be-46c4-b67e-5119c8c5de53/submissions/231ea80e-bef3-4e96-90a4-07efc80eb523/minimap2/20583855-420b-4027-97a0-b521fa94934e/call-alignAndSortBAM/HG00673.bam";

    Timer t;

    while (true){
        authenticator.update();

        for_read_in_bam_region(bam_path, region_string, [&](Sequence& s){
            cerr << t << ' ' << s.name << '\n';
            // Stop after the first read
            return;
        });

        sleep(30);
    }

}
