//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: CASPT2_densityq3.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <bagel_config.h>
#ifdef COMPILE_SMITH


#include <src/smith/caspt2/CASPT2.h>
#include <src/smith/caspt2/CASPT2_tasks9.h>
#include <src/smith/caspt2/CASPT2_tasks10.h>
#include <src/smith/caspt2/CASPT2_tasks11.h>

using namespace std;
using namespace bagel;
using namespace bagel::SMITH;

void  CASPT2::CASPT2::make_densityq3(shared_ptr<Queue> densityq, shared_ptr<Task> task284, const bool diagonal) {
  array<shared_ptr<const IndexRange>,3> pindex = {{rclosed_, ractive_, rvirt_}};

  vector<IndexRange> I500_index = {active_, virt_};
  auto I500 = make_shared<Tensor>(I500_index);
  auto tensor424 = vector<shared_ptr<Tensor>>{den2, I500};
  auto task424 = make_shared<Task424>(tensor424, pindex);
  task424->add_dep(task284);
  densityq->add_task(task424);

  vector<IndexRange> I501_index = {closed_, active_, active_, active_};
  auto I501 = make_shared<Tensor>(I501_index);
  auto tensor425 = vector<shared_ptr<Tensor>>{I500, t2, I501};
  auto task425 = make_shared<Task425>(tensor425, pindex);
  task424->add_dep(task425);
  task425->add_dep(task284);
  densityq->add_task(task425);

  auto tensor426 = vector<shared_ptr<Tensor>>{I501, Gamma6_(), t2};
  auto task426 = make_shared<Task426>(tensor426, pindex);
  task425->add_dep(task426);
  task426->add_dep(task284);
  densityq->add_task(task426);

  vector<IndexRange> I653_index = {virt_, active_, active_, active_};
  auto I653 = make_shared<Tensor>(I653_index);
  auto tensor427 = vector<shared_ptr<Tensor>>{I500, t2, I653};
  auto task427 = make_shared<Task427>(tensor427, pindex);
  task424->add_dep(task427);
  task427->add_dep(task284);
  densityq->add_task(task427);

  auto tensor428 = vector<shared_ptr<Tensor>>{I653, Gamma59_(), t2};
  auto task428 = make_shared<Task428>(tensor428, pindex);
  task427->add_dep(task428);
  task428->add_dep(task284);
  densityq->add_task(task428);

  vector<IndexRange> I512_index = {closed_, closed_};
  auto I512 = make_shared<Tensor>(I512_index);
  auto tensor429 = vector<shared_ptr<Tensor>>{den2, I512};
  auto task429 = make_shared<Task429>(tensor429, pindex);
  task429->add_dep(task284);
  densityq->add_task(task429);

  vector<IndexRange> I513_index = {virt_, closed_, active_, active_};
  auto I513 = make_shared<Tensor>(I513_index);
  auto tensor430 = vector<shared_ptr<Tensor>>{I512, t2, I513};
  auto task430 = make_shared<Task430>(tensor430, pindex);
  task429->add_dep(task430);
  task430->add_dep(task284);
  densityq->add_task(task430);

  auto tensor431 = vector<shared_ptr<Tensor>>{I513, Gamma35_(), t2};
  auto task431 = make_shared<Task431>(tensor431, pindex);
  task430->add_dep(task431);
  task431->add_dep(task284);
  densityq->add_task(task431);

  vector<IndexRange> I522_index = {virt_, closed_, active_, active_};
  auto I522 = make_shared<Tensor>(I522_index);
  auto tensor432 = vector<shared_ptr<Tensor>>{I512, t2, I522};
  auto task432 = make_shared<Task432>(tensor432, pindex);
  task429->add_dep(task432);
  task432->add_dep(task284);
  densityq->add_task(task432);

  auto tensor433 = vector<shared_ptr<Tensor>>{I522, Gamma35_(), t2};
  auto task433 = make_shared<Task433>(tensor433, pindex);
  task432->add_dep(task433);
  task433->add_dep(task284);
  densityq->add_task(task433);

  vector<IndexRange> I515_index = {virt_, virt_};
  auto I515 = make_shared<Tensor>(I515_index);
  auto tensor434 = vector<shared_ptr<Tensor>>{den2, I515};
  auto task434 = make_shared<Task434>(tensor434, pindex);
  task434->add_dep(task284);
  densityq->add_task(task434);

  vector<IndexRange> I516_index = {virt_, closed_, active_, active_};
  auto I516 = make_shared<Tensor>(I516_index);
  auto tensor435 = vector<shared_ptr<Tensor>>{I515, t2, I516};
  auto task435 = make_shared<Task435>(tensor435, pindex);
  task434->add_dep(task435);
  task435->add_dep(task284);
  densityq->add_task(task435);

  auto tensor436 = vector<shared_ptr<Tensor>>{I516, Gamma35_(), t2};
  auto task436 = make_shared<Task436>(tensor436, pindex);
  task435->add_dep(task436);
  task436->add_dep(task284);
  densityq->add_task(task436);

  vector<IndexRange> I525_index = {virt_, closed_, active_, active_};
  auto I525 = make_shared<Tensor>(I525_index);
  auto tensor437 = vector<shared_ptr<Tensor>>{I515, t2, I525};
  auto task437 = make_shared<Task437>(tensor437, pindex);
  task434->add_dep(task437);
  task437->add_dep(task284);
  densityq->add_task(task437);

  auto tensor438 = vector<shared_ptr<Tensor>>{I525, Gamma35_(), t2};
  auto task438 = make_shared<Task438>(tensor438, pindex);
  task437->add_dep(task438);
  task438->add_dep(task284);
  densityq->add_task(task438);

  vector<IndexRange> I662_index = {virt_, virt_, active_, active_};
  auto I662 = make_shared<Tensor>(I662_index);
  auto tensor439 = vector<shared_ptr<Tensor>>{I515, t2, I662};
  auto task439 = make_shared<Task439>(tensor439, pindex);
  task434->add_dep(task439);
  task439->add_dep(task284);
  densityq->add_task(task439);

  auto tensor440 = vector<shared_ptr<Tensor>>{I662, Gamma60_(), t2};
  auto task440 = make_shared<Task440>(tensor440, pindex);
  task439->add_dep(task440);
  task440->add_dep(task284);
  densityq->add_task(task440);

  vector<IndexRange> I527_index = {active_, closed_};
  auto I527 = make_shared<Tensor>(I527_index);
  auto tensor441 = vector<shared_ptr<Tensor>>{den2, I527};
  auto task441 = make_shared<Task441>(tensor441, pindex);
  task441->add_dep(task284);
  densityq->add_task(task441);

  vector<IndexRange> I528_index = {virt_, active_, active_, active_};
  auto I528 = make_shared<Tensor>(I528_index);
  auto tensor442 = vector<shared_ptr<Tensor>>{I527, t2, I528};
  auto task442 = make_shared<Task442>(tensor442, pindex);
  task441->add_dep(task442);
  task442->add_dep(task284);
  densityq->add_task(task442);

  auto tensor443 = vector<shared_ptr<Tensor>>{I528, Gamma51_(), t2};
  auto task443 = make_shared<Task443>(tensor443, pindex);
  task442->add_dep(task443);
  task443->add_dep(task284);
  densityq->add_task(task443);

  vector<IndexRange> I551_index = {virt_, virt_};
  auto I551 = make_shared<Tensor>(I551_index);
  auto tensor444 = vector<shared_ptr<Tensor>>{den2, I551};
  auto task444 = make_shared<Task444>(tensor444, pindex);
  task444->add_dep(task284);
  densityq->add_task(task444);

  vector<IndexRange> I552_index = {virt_, active_, active_, active_};
  auto I552 = make_shared<Tensor>(I552_index);
  auto tensor445 = vector<shared_ptr<Tensor>>{I551, t2, I552};
  auto task445 = make_shared<Task445>(tensor445, pindex);
  task444->add_dep(task445);
  task445->add_dep(task284);
  densityq->add_task(task445);

  auto tensor446 = vector<shared_ptr<Tensor>>{I552, Gamma59_(), t2};
  auto task446 = make_shared<Task446>(tensor446, pindex);
  task445->add_dep(task446);
  task446->add_dep(task284);
  densityq->add_task(task446);

  vector<IndexRange> I563_index = {virt_, active_};
  auto I563 = make_shared<Tensor>(I563_index);
  auto tensor447 = vector<shared_ptr<Tensor>>{den2, I563};
  auto task447 = make_shared<Task447>(tensor447, pindex);
  task447->add_dep(task284);
  densityq->add_task(task447);

  vector<IndexRange> I564_index = {virt_, active_};
  auto I564 = make_shared<Tensor>(I564_index);
  auto tensor448 = vector<shared_ptr<Tensor>>{I563, Gamma16_(), I564};
  auto task448 = make_shared<Task448>(tensor448, pindex);
  task447->add_dep(task448);
  task448->add_dep(task284);
  densityq->add_task(task448);

  auto tensor449 = vector<shared_ptr<Tensor>>{I564, t2};
  auto task449 = make_shared<Task449>(tensor449, pindex);
  task448->add_dep(task449);
  task449->add_dep(task284);
  densityq->add_task(task449);

  vector<IndexRange> I566_index = {virt_, active_};
  auto I566 = make_shared<Tensor>(I566_index);
  auto tensor450 = vector<shared_ptr<Tensor>>{den2, I566};
  auto task450 = make_shared<Task450>(tensor450, pindex);
  task450->add_dep(task284);
  densityq->add_task(task450);

  vector<IndexRange> I567_index = {virt_, active_};
  auto I567 = make_shared<Tensor>(I567_index);
  auto tensor451 = vector<shared_ptr<Tensor>>{I566, Gamma16_(), I567};
  auto task451 = make_shared<Task451>(tensor451, pindex);
  task450->add_dep(task451);
  task451->add_dep(task284);
  densityq->add_task(task451);

  auto tensor452 = vector<shared_ptr<Tensor>>{I567, t2};
  auto task452 = make_shared<Task452>(tensor452, pindex);
  task451->add_dep(task452);
  task452->add_dep(task284);
  densityq->add_task(task452);

  vector<IndexRange> I572_index = {virt_, closed_};
  auto I572 = make_shared<Tensor>(I572_index);
  auto tensor453 = vector<shared_ptr<Tensor>>{den2, I572};
  auto task453 = make_shared<Task453>(tensor453, pindex);
  task453->add_dep(task284);
  densityq->add_task(task453);

  vector<IndexRange> I573_index = {virt_, closed_};
  auto I573 = make_shared<Tensor>(I573_index);
  auto tensor454 = vector<shared_ptr<Tensor>>{I572, t2, I573};
  auto task454 = make_shared<Task454>(tensor454, pindex);
  task453->add_dep(task454);
  task454->add_dep(task284);
  densityq->add_task(task454);

  vector<IndexRange> I574_index = {active_, virt_, closed_, active_};
  auto I574 = make_shared<Tensor>(I574_index);
  auto tensor455 = vector<shared_ptr<Tensor>>{I573, Gamma38_(), I574};
  auto task455 = make_shared<Task455>(tensor455, pindex);
  task454->add_dep(task455);
  task455->add_dep(task284);
  densityq->add_task(task455);

  auto tensor456 = vector<shared_ptr<Tensor>>{I574, t2};
  auto task456 = make_shared<Task456>(tensor456, pindex);
  task455->add_dep(task456);
  task456->add_dep(task284);
  densityq->add_task(task456);

  vector<IndexRange> I581_index = {active_, active_};
  auto I581 = make_shared<Tensor>(I581_index);
  auto tensor457 = vector<shared_ptr<Tensor>>{den2, I581};
  auto task457 = make_shared<Task457>(tensor457, pindex);
  task457->add_dep(task284);
  densityq->add_task(task457);

  vector<IndexRange> I582_index;
  auto I582 = make_shared<Tensor>(I582_index);
  auto tensor458 = vector<shared_ptr<Tensor>>{I581, Gamma38_(), I582};
  auto task458 = make_shared<Task458>(tensor458, pindex);
  task457->add_dep(task458);
  task458->add_dep(task284);
  densityq->add_task(task458);

  auto tensor459 = vector<shared_ptr<Tensor>>{I582, t2};
  auto task459 = make_shared<Task459>(tensor459, pindex);
  task458->add_dep(task459);
  task459->add_dep(task284);
  densityq->add_task(task459);

  auto tensor460 = vector<shared_ptr<Tensor>>{I582, t2};
  auto task460 = make_shared<Task460>(tensor460, pindex);
  task458->add_dep(task460);
  task460->add_dep(task284);
  densityq->add_task(task460);

  shared_ptr<Tensor> I587;
  if (diagonal) {
    vector<IndexRange> I587_index = {closed_, closed_};
    I587 = make_shared<Tensor>(I587_index);
  }
  shared_ptr<Task461> task461;
  if (diagonal) {
    auto tensor461 = vector<shared_ptr<Tensor>>{den2, I587};
    task461 = make_shared<Task461>(tensor461, pindex);
    task461->add_dep(task284);
    densityq->add_task(task461);
  }

  shared_ptr<Task462> task462;
  if (diagonal) {
    auto tensor462 = vector<shared_ptr<Tensor>>{I587, t2};
    task462 = make_shared<Task462>(tensor462, pindex);
    task461->add_dep(task462);
    task462->add_dep(task284);
    densityq->add_task(task462);
  }

  shared_ptr<Task463> task463;
  if (diagonal) {
    auto tensor463 = vector<shared_ptr<Tensor>>{I587, t2};
    task463 = make_shared<Task463>(tensor463, pindex);
    task461->add_dep(task463);
    task463->add_dep(task284);
    densityq->add_task(task463);
  }

  shared_ptr<Tensor> I591;
  if (diagonal) {
    vector<IndexRange> I591_index = {virt_, virt_};
    I591 = make_shared<Tensor>(I591_index);
  }
  shared_ptr<Task464> task464;
  if (diagonal) {
    auto tensor464 = vector<shared_ptr<Tensor>>{den2, I591};
    task464 = make_shared<Task464>(tensor464, pindex);
    task464->add_dep(task284);
    densityq->add_task(task464);
  }

  shared_ptr<Task465> task465;
  if (diagonal) {
    auto tensor465 = vector<shared_ptr<Tensor>>{I591, t2};
    task465 = make_shared<Task465>(tensor465, pindex);
    task464->add_dep(task465);
    task465->add_dep(task284);
    densityq->add_task(task465);
  }

  shared_ptr<Task466> task466;
  if (diagonal) {
    auto tensor466 = vector<shared_ptr<Tensor>>{I591, t2};
    task466 = make_shared<Task466>(tensor466, pindex);
    task464->add_dep(task466);
    task466->add_dep(task284);
    densityq->add_task(task466);
  }

  vector<IndexRange> I595_index = {closed_, active_};
  auto I595 = make_shared<Tensor>(I595_index);
  auto tensor467 = vector<shared_ptr<Tensor>>{den2, I595};
  auto task467 = make_shared<Task467>(tensor467, pindex);
  task467->add_dep(task284);
  densityq->add_task(task467);

  vector<IndexRange> I596_index = {closed_, active_};
  auto I596 = make_shared<Tensor>(I596_index);
  auto tensor468 = vector<shared_ptr<Tensor>>{I595, Gamma38_(), I596};
  auto task468 = make_shared<Task468>(tensor468, pindex);
  task467->add_dep(task468);
  task468->add_dep(task284);
  densityq->add_task(task468);

  auto tensor469 = vector<shared_ptr<Tensor>>{I596, t2};
  auto task469 = make_shared<Task469>(tensor469, pindex);
  task468->add_dep(task469);
  task469->add_dep(task284);
  densityq->add_task(task469);

  auto tensor470 = vector<shared_ptr<Tensor>>{I596, t2};
  auto task470 = make_shared<Task470>(tensor470, pindex);
  task468->add_dep(task470);
  task470->add_dep(task284);
  densityq->add_task(task470);

  vector<IndexRange> I601_index = {active_, virt_};
  auto I601 = make_shared<Tensor>(I601_index);
  auto tensor471 = vector<shared_ptr<Tensor>>{den2, I601};
  auto task471 = make_shared<Task471>(tensor471, pindex);
  task471->add_dep(task284);
  densityq->add_task(task471);

  vector<IndexRange> I602_index = {virt_, closed_, active_, active_};
  auto I602 = make_shared<Tensor>(I602_index);
  auto tensor472 = vector<shared_ptr<Tensor>>{I601, t2, I602};
  auto task472 = make_shared<Task472>(tensor472, pindex);
  task471->add_dep(task472);
  task472->add_dep(task284);
  densityq->add_task(task472);

  vector<IndexRange> I603_index = {active_, virt_, closed_, active_};
  auto I603 = make_shared<Tensor>(I603_index);
  auto tensor473 = vector<shared_ptr<Tensor>>{I602, Gamma35_(), I603};
  auto task473 = make_shared<Task473>(tensor473, pindex);
  task472->add_dep(task473);
  task473->add_dep(task284);
  densityq->add_task(task473);

  auto tensor474 = vector<shared_ptr<Tensor>>{I603, t2};
  auto task474 = make_shared<Task474>(tensor474, pindex);
  task473->add_dep(task474);
  task474->add_dep(task284);
  densityq->add_task(task474);

  vector<IndexRange> I604_index = {active_, virt_};
  auto I604 = make_shared<Tensor>(I604_index);
  auto tensor475 = vector<shared_ptr<Tensor>>{den2, I604};
  auto task475 = make_shared<Task475>(tensor475, pindex);
  task475->add_dep(task284);
  densityq->add_task(task475);

  vector<IndexRange> I605_index = {virt_, closed_, active_, active_};
  auto I605 = make_shared<Tensor>(I605_index);
  auto tensor476 = vector<shared_ptr<Tensor>>{I604, t2, I605};
  auto task476 = make_shared<Task476>(tensor476, pindex);
  task475->add_dep(task476);
  task476->add_dep(task284);
  densityq->add_task(task476);

  auto tensor477 = vector<shared_ptr<Tensor>>{I605, Gamma32_(), t2};
  auto task477 = make_shared<Task477>(tensor477, pindex);
  task476->add_dep(task477);
  task477->add_dep(task284);
  densityq->add_task(task477);

  auto tensor478 = vector<shared_ptr<Tensor>>{I605, Gamma35_(), t2};
  auto task478 = make_shared<Task478>(tensor478, pindex);
  task476->add_dep(task478);
  task478->add_dep(task284);
  densityq->add_task(task478);

  vector<IndexRange> I613_index = {closed_, virt_};
  auto I613 = make_shared<Tensor>(I613_index);
  auto tensor479 = vector<shared_ptr<Tensor>>{den2, I613};
  auto task479 = make_shared<Task479>(tensor479, pindex);
  task479->add_dep(task284);
  densityq->add_task(task479);

  vector<IndexRange> I614_index = {virt_, active_};
  auto I614 = make_shared<Tensor>(I614_index);
  auto tensor480 = vector<shared_ptr<Tensor>>{I613, t2, I614};
  auto task480 = make_shared<Task480>(tensor480, pindex);
  task479->add_dep(task480);
  task480->add_dep(task284);
  densityq->add_task(task480);

  auto tensor481 = vector<shared_ptr<Tensor>>{I614, Gamma60_(), t2};
  auto task481 = make_shared<Task481>(tensor481, pindex);
  task480->add_dep(task481);
  task481->add_dep(task284);
  densityq->add_task(task481);

  vector<IndexRange> I616_index = {virt_, closed_};
  auto I616 = make_shared<Tensor>(I616_index);
  auto tensor482 = vector<shared_ptr<Tensor>>{den2, I616};
  auto task482 = make_shared<Task482>(tensor482, pindex);
  task482->add_dep(task284);
  densityq->add_task(task482);

  vector<IndexRange> I617_index = {virt_, active_};
  auto I617 = make_shared<Tensor>(I617_index);
  auto tensor483 = vector<shared_ptr<Tensor>>{I616, t2, I617};
  auto task483 = make_shared<Task483>(tensor483, pindex);
  task482->add_dep(task483);
  task483->add_dep(task284);
  densityq->add_task(task483);

  auto tensor484 = vector<shared_ptr<Tensor>>{I617, Gamma60_(), t2};
  auto task484 = make_shared<Task484>(tensor484, pindex);
  task483->add_dep(task484);
  task484->add_dep(task284);
  densityq->add_task(task484);

  vector<IndexRange> I619_index = {closed_, active_};
  auto I619 = make_shared<Tensor>(I619_index);
  auto tensor485 = vector<shared_ptr<Tensor>>{den2, I619};
  auto task485 = make_shared<Task485>(tensor485, pindex);
  task485->add_dep(task284);
  densityq->add_task(task485);

  vector<IndexRange> I620_index = {active_, closed_};
  auto I620 = make_shared<Tensor>(I620_index);
  auto tensor486 = vector<shared_ptr<Tensor>>{I619, Gamma38_(), I620};
  auto task486 = make_shared<Task486>(tensor486, pindex);
  task485->add_dep(task486);
  task486->add_dep(task284);
  densityq->add_task(task486);

  auto tensor487 = vector<shared_ptr<Tensor>>{I620, t2};
  auto task487 = make_shared<Task487>(tensor487, pindex);
  task486->add_dep(task487);
  task487->add_dep(task284);
  densityq->add_task(task487);

  auto tensor488 = vector<shared_ptr<Tensor>>{I620, t2};
  auto task488 = make_shared<Task488>(tensor488, pindex);
  task486->add_dep(task488);
  task488->add_dep(task284);
  densityq->add_task(task488);

  vector<IndexRange> I631_index = {closed_, closed_};
  auto I631 = make_shared<Tensor>(I631_index);
  auto tensor489 = vector<shared_ptr<Tensor>>{den2, I631};
  auto task489 = make_shared<Task489>(tensor489, pindex);
  task489->add_dep(task284);
  densityq->add_task(task489);

  vector<IndexRange> I632_index = {virt_, closed_, virt_, active_};
  auto I632 = make_shared<Tensor>(I632_index);
  auto tensor490 = vector<shared_ptr<Tensor>>{I631, t2, I632};
  auto task490 = make_shared<Task490>(tensor490, pindex);
  task489->add_dep(task490);
  task490->add_dep(task284);
  densityq->add_task(task490);

  auto tensor491 = vector<shared_ptr<Tensor>>{I632, Gamma38_(), t2};
  auto task491 = make_shared<Task491>(tensor491, pindex);
  task490->add_dep(task491);
  task491->add_dep(task284);
  densityq->add_task(task491);

  vector<IndexRange> I635_index = {virt_, closed_, virt_, active_};
  auto I635 = make_shared<Tensor>(I635_index);
  auto tensor492 = vector<shared_ptr<Tensor>>{I631, t2, I635};
  auto task492 = make_shared<Task492>(tensor492, pindex);
  task489->add_dep(task492);
  task492->add_dep(task284);
  densityq->add_task(task492);

  auto tensor493 = vector<shared_ptr<Tensor>>{I635, Gamma38_(), t2};
  auto task493 = make_shared<Task493>(tensor493, pindex);
  task492->add_dep(task493);
  task493->add_dep(task284);
  densityq->add_task(task493);

  vector<IndexRange> I637_index = {virt_, virt_};
  auto I637 = make_shared<Tensor>(I637_index);
  auto tensor494 = vector<shared_ptr<Tensor>>{den2, I637};
  auto task494 = make_shared<Task494>(tensor494, pindex);
  task494->add_dep(task284);
  densityq->add_task(task494);

  vector<IndexRange> I638_index = {virt_, closed_, virt_, active_};
  auto I638 = make_shared<Tensor>(I638_index);
  auto tensor495 = vector<shared_ptr<Tensor>>{I637, t2, I638};
  auto task495 = make_shared<Task495>(tensor495, pindex);
  task494->add_dep(task495);
  task495->add_dep(task284);
  densityq->add_task(task495);

  auto tensor496 = vector<shared_ptr<Tensor>>{I638, Gamma38_(), t2};
  auto task496 = make_shared<Task496>(tensor496, pindex);
  task495->add_dep(task496);
  task496->add_dep(task284);
  densityq->add_task(task496);

  vector<IndexRange> I641_index = {virt_, closed_, virt_, active_};
  auto I641 = make_shared<Tensor>(I641_index);
  auto tensor497 = vector<shared_ptr<Tensor>>{I637, t2, I641};
  auto task497 = make_shared<Task497>(tensor497, pindex);
  task494->add_dep(task497);
  task497->add_dep(task284);
  densityq->add_task(task497);

  auto tensor498 = vector<shared_ptr<Tensor>>{I641, Gamma38_(), t2};
  auto task498 = make_shared<Task498>(tensor498, pindex);
  task497->add_dep(task498);
  task498->add_dep(task284);
  densityq->add_task(task498);

  vector<IndexRange> I643_index = {virt_, virt_};
  auto I643 = make_shared<Tensor>(I643_index);
  auto tensor499 = vector<shared_ptr<Tensor>>{den2, I643};
  auto task499 = make_shared<Task499>(tensor499, pindex);
  task499->add_dep(task284);
  densityq->add_task(task499);

  vector<IndexRange> I644_index = {virt_, closed_, virt_, active_};
  auto I644 = make_shared<Tensor>(I644_index);
  auto tensor500 = vector<shared_ptr<Tensor>>{I643, t2, I644};
  auto task500 = make_shared<Task500>(tensor500, pindex);
  task499->add_dep(task500);
  task500->add_dep(task284);
  densityq->add_task(task500);

  auto tensor501 = vector<shared_ptr<Tensor>>{I644, Gamma38_(), t2};
  auto task501 = make_shared<Task501>(tensor501, pindex);
  task500->add_dep(task501);
  task501->add_dep(task284);
  densityq->add_task(task501);

  vector<IndexRange> I647_index = {virt_, closed_, virt_, active_};
  auto I647 = make_shared<Tensor>(I647_index);
  auto tensor502 = vector<shared_ptr<Tensor>>{I643, t2, I647};
  auto task502 = make_shared<Task502>(tensor502, pindex);
  task499->add_dep(task502);
  task502->add_dep(task284);
  densityq->add_task(task502);

  auto tensor503 = vector<shared_ptr<Tensor>>{I647, Gamma38_(), t2};
  auto task503 = make_shared<Task503>(tensor503, pindex);
  task502->add_dep(task503);
  task503->add_dep(task284);
  densityq->add_task(task503);

  vector<IndexRange> I649_index = {active_, closed_};
  auto I649 = make_shared<Tensor>(I649_index);
  auto tensor504 = vector<shared_ptr<Tensor>>{den2, I649};
  auto task504 = make_shared<Task504>(tensor504, pindex);
  task504->add_dep(task284);
  densityq->add_task(task504);

  vector<IndexRange> I650_index = {virt_, virt_, active_, active_};
  auto I650 = make_shared<Tensor>(I650_index);
  auto tensor505 = vector<shared_ptr<Tensor>>{I649, t2, I650};
  auto task505 = make_shared<Task505>(tensor505, pindex);
  task504->add_dep(task505);
  task505->add_dep(task284);
  densityq->add_task(task505);

  auto tensor506 = vector<shared_ptr<Tensor>>{I650, Gamma60_(), t2};
  auto task506 = make_shared<Task506>(tensor506, pindex);
  task505->add_dep(task506);
  task506->add_dep(task284);
  densityq->add_task(task506);
}

#endif
