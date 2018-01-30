import io;
import sys;
import files;

string emews_root = getenv("EMEWS_PROJECT_ROOT");
string turbine_output = getenv("TURBINE_OUTPUT");

app (file out, file err) run_model (file shfile, string param_line, string instance)
{
    "bash" shfile param_line emews_root instance @stdout=out @stderr=err;
}

// call this to create any required directories
app (void o) make_dir(string dirname) {
  "mkdir" "-p" dirname;
}

// anything that need to be done prior to a model runs
// (e.g. file creation) can be done here
//app (void o) run_prerequisites() {
//
//}

//run_prerequisites() => {
  file model_sh = input(emews_root+"/scripts/cancer-emews.sh");
  file upf = input(argv("f"));
  string upf_lines[] = file_lines(upf);
  foreach s,i in upf_lines {
    string instance = "%s/instance_%i/" % (turbine_output, i+1);
    make_dir(instance) => {
      file out <instance+"out.txt">;
      file err <instance+"err.txt">;
      (out,err) = run_model(model_sh, s, instance);
    }
  }
//}
