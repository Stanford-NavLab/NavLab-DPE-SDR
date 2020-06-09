#!/bin/sh

# run sudo su before hand to be in root mode

echo 0 > /sys/kernel/debug/clock/override.emc/state
echo "EMC override state: `cat /sys/kernel/debug/clock/override.emc/state`"
echo "EMC rate: `cat /sys/kernel/debug/clock/emc/rate`"

echo 0 > /sys/kernel/debug/clock/override.gbus/state
echo "GPU override state: `cat /sys/kernel/debug/clock/override.gbus/state`"
echo "GPU rate: `cat /sys/kernel/debug/clock/gbus/rate`"

echo interactive > /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor
echo "CPU Scaling governor: `cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor`"
for i in 0 1 2 3 ; do
	echo "CPU${i}: `cat /sys/devices/system/cpu/cpu${i}/cpufreq/scaling_cur_freq`"
done

# turn off fan
echo "Disabling fan..."
if [ ! -w /sys/kernel/debug/tegra_fan/target_pwm ] ; then
	echo "Cannot set fan -- exiting..."
fi
echo 0 > /sys/kernel/debug/tegra_fan/target_pwm

