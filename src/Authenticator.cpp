#include "Authenticator.hpp"
#include "curl/curl.h"
#include <functional>
#include <exception>
#include <thread>

using std::function;
using std::exception;
using std::this_thread::sleep_for;

namespace sv_merge{



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
    }
}


void parse_json_line(const string& line, pair<string,string>& result){
    result.first.clear();
    result.second.clear();
    size_t n_quotes = 0;
    size_t n_colons = 0;

    for (auto c: line){
        if (c == '\"'){
            n_quotes++;
            continue;
        }
        if (c == ':'){
            if (n_quotes < 2 or n_quotes > 3){
                throw runtime_error("ERROR: attempt to parse nested value in JSON line: " + line);
            }

            n_colons++;
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


int64_t get_remaining_seconds(const string& token){
    string url = "https://oauth2.googleapis.com/tokeninfo?access_token=" + token;
    string result;
    curl_get(url, result);

    int64_t remaining = -1;
    for_item_in_json(result, [&](const pair<string,string>& item){
        if (item.second == "invalid_token"){
            remaining = -1;
            cerr << "TOKEN EXPIRED" << '\n';
            return;
        }
        else if (item.first == "error"){
            throw runtime_error("ERROR: cannot acquire token info: " + item.second);
        }
        else if (item.first == "expires_in"){
            remaining = stoll(item.second);
            cerr << "SUCCESS: time remaining " << remaining << " sec" << '\n';
            return;
        }
    });

    return remaining;
}


void GoogleAuthenticator::update() {
    m.lock();

    // Don't bother querying the info server if we aren't withing 60 seconds of the expiration
    if (get_current_time() < expiration - seconds(60)){
        m.unlock();
        return;
    }

    // Request info from the server and update the expiration timepoint
    auto s = get_remaining_seconds(token);
    expiration = get_current_time() + seconds(s);

    int64_t delay_sec = 0;
    int64_t n_retries = 0;

    // Attempt to update the token, with some number of retries, and an exponential backoff period
    while(s < 0){
        cerr << "Waiting " << delay_sec << " sec" << '\n';
        sleep_for(seconds(delay_sec));

        cerr << "Updating token" << '\n';
        string command = "gcloud auth print-access-token";

        auto t = get_current_time();

        // Will throw error if fails, otherwise returns token
        run_command(command, token);

        string name = "GCS_OAUTH_TOKEN";
        string env;

        // Update the environment with the token
        auto error_code = setenv("GCS_OAUTH_TOKEN", token.c_str(), 1);

        auto result = getenv(name.data());

        if (result == nullptr or error_code != 0) {
            cerr << "error_code: " << error_code << '\n';
            throw runtime_error("ERROR: environment variable not set");
        } else {
            env = result;
        }

        n_retries++;

        if (n_retries >= 5){
            throw runtime_error("ERROR: authentication failed after " + to_string(n_retries) + " retries");
        }

        s = get_remaining_seconds(token);
        expiration = get_current_time() + seconds(s);

        if (delay_sec == 0){
            delay_sec = 1;
        }
        else{
            delay_sec *= 2;
        }
    }

    m.unlock();
}


void GoogleAuthenticator::try_with_authentication(int64_t n_retries, const function<void()>& f){
    update();
    int64_t n = 0;
    bool success = false;

    while (not success and n < n_retries) {
        try {
            f();
        }
        catch (exception& e) {
            cerr << e.what() << '\n';
            cerr << "Authenticating..." << '\n';
            update();
            cerr << "Retrying..." << '\n';
            n++;
        }

        success = true;
    }
}


}
