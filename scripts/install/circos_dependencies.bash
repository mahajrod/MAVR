#!/usr/bin/env bash

sudo apt-get install -y cpanminus libgd2-dev

sudo cpanm Config::General GD Font::TTF::Font GD::Polyline Math::Bezier Params::Validate Math::VecStat Readonly Statistics::Basic Regexp::Common SVG Set::IntSpan Text::Format
