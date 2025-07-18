#!/bin/bash

urls=(
  "https://www.dropbox.com/scl/fi/az5665h18b9ancc8fq0c4/chi_exact_distributions.7z.part_000?rlkey=hjii4y6vzk99txv21olzl0kki&st=89vi3v9f&dl=1"
  "https://www.dropbox.com/scl/fi/eih9tg0ygk1897uotw4r8/chi_exact_distributions.7z.part_001?rlkey=hcaho1m729tk3ynd2l14rm4ye&st=d95gn7w4&dl=1"
  "https://www.dropbox.com/scl/fi/e7iwn1s6v53qcpl5mttcz/chi_exact_distributions.7z.part_002?rlkey=xrfrs71v8l27swqzmz9zu7xrj&st=1x3gbcjg&dl=1"
  "https://www.dropbox.com/scl/fi/700zb3w49firypn7eeidt/chi_exact_distributions.7z.part_003?rlkey=5szx5ga9x8flqc2vqdcof5n2x&st=dcwwpxif&dl=1"
  "https://www.dropbox.com/scl/fi/h786v1h6bplkdq19huzy9/chi_exact_distributions.7z.part_004?rlkey=emu4q1djbithqjhz6uw0g8c8k&st=60zpbtpu&dl=1"
  "https://www.dropbox.com/scl/fi/eq8ixxqwswfc4n1f6jtmn/chi_exact_distributions.7z.part_005?rlkey=2r58k9ccz7rbcn5a44t2oafg8&st=owwnm605&dl=1"
  "https://www.dropbox.com/scl/fi/13qrqj2mndf7qgbit6t9a/chi_exact_distributions.7z.part_006?rlkey=98iwlo4wq8l7pd96u71yjjzuw&st=efw4pqxi&dl=1"
  "https://www.dropbox.com/scl/fi/zveku7hwvyblpiwdxejno/chi_exact_distributions.7z.part_007?rlkey=ilu4aofvoj3y4wjm16lhd6tpb&st=wg31qnfz&dl=1"
)

output="chi_exact_distributions.7z"

rm -f "$output"

for url in "${urls[@]}"; do
    echo "Downloading and appending from $url..."
    wget -q -O - "$url" >> "$output"
done

echo "All parts downloaded and appended into $output."
