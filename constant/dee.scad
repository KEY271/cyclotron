translate([0, 0, 1]) {
    difference() {
        difference() {
            cylinder(h=15, r=25, $fn=100);
            translate([25, 0, 0]) {
                cube(50, center=true);
            }
        }
        translate([0, 0, 3]) {
            cylinder(h=9, r=24, $fn=100);
        }
    }
}
