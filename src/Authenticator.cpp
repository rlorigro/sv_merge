#include "Authenticator.hpp"


namespace sv_merge{

void GoogleAuthenticator::update() {
    // Check if expires within 30sec
    if (expiration - seconds(30) < get_current_time()) {
        cerr << "Updating token" << '\n';
        string command = "gcloud auth print-access-token";
        string token;

        // Will throw error if fails, otherwise returns token
        auto t = get_current_time();
        run_command(command, token);
        expiration = t + token_lifetime;

        string name = "GCS_OAUTH_TOKEN";
        string env;

        // Update the environment with the token
        auto error_code = setenv("GCS_OAUTH_TOKEN", token.c_str(), 1);
        cerr << "result: " << error_code << '\n';

        auto result = getenv(name.data());

        if (result == nullptr or error_code != 0) {
            throw runtime_error("ERROR: environment variable not set");
        } else {
            env = result;
        }
    }
}

}
